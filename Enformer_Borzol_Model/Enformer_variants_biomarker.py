"""
Enformer Multi-Sample Analysis for Biomarker Extraction

This script demonstrates how to:
1. Process multiple sample variants from a VCF file
2. Compare each sample to a reference genome using Enformer
3. Identify consistent regulatory changes across samples
4. Extract potential biomarkers from regulatory patterns
5. Cluster samples based on their regulatory profiles
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tensorflow as tf
import tensorflow_hub as hub
from kipoiseq import Interval
from kipoiseq.extractors import FastaStringExtractor
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
import scipy.stats as stats
import vcf  # PyVCF package

# Constants
SEQUENCE_LENGTH = 393_216  # Input sequence length for Enformer
BIN_SIZE = 128  # Size of bins for prediction
CENTRAL_BINS = 896  # Number of central bins for prediction

class EnformerMultiSampleAnalyzer:
    """Class for analyzing multiple samples using Enformer."""
    
    def __init__(self, model_path='https://tfhub.dev/deepmind/enformer/1', 
                 genome_fasta_path='hg38.fa'):
        """Initialize the analyzer.
        
        Args:
            model_path: Path to the Enformer model.
            genome_fasta_path: Path to the reference genome FASTA file.
        """
        print("Loading Enformer model...")
        self.model = hub.load(model_path)
        print("Model loaded successfully!")
        
        # Check if genome file exists before loading
        if os.path.exists(genome_fasta_path):
            print(f"Loading reference genome from {genome_fasta_path}...")
            self.fasta_extractor = FastaStringExtractor(genome_fasta_path)
        else:
            print(f"Warning: Genome file {genome_fasta_path} not found.")
            self.fasta_extractor = None
            
        # Load target information
        self.targets_df = self._load_target_info()
        
    def _load_target_info(self):
        """Load information about the target tracks."""
        try:
            targets_url = 'https://raw.githubusercontent.com/calico/basenji/0.5/manuscripts/cross2020/targets_human.txt'
            df = pd.read_csv(targets_url, sep='\t')
            print(f"Loaded {len(df)} target tracks")
            return df
        except:
            print("Could not load targets info. Using mock data.")
            # Create mock data
            return pd.DataFrame({
                'index': range(5313),
                'description': [f'Target {i}' for i in range(5313)],
                'file': ['file.txt'] * 5313,
                'type': ['CAGE'] * 1000 + ['DNASE'] * 1000 + ['CHIP'] * 3313,
                'cell_type': ['cell_type'] * 5313
            })
    
    def one_hot_encode(self, sequence):
        """Convert DNA sequence to one-hot encoding.
        
        Args:
            sequence: DNA sequence as a string.
            
        Returns:
            One-hot encoded sequence as numpy array.
        """
        sequence = sequence.upper()
        nucleotide_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        one_hot = np.zeros((len(sequence), 4), dtype=np.float32)
        
        for i, nucleotide in enumerate(sequence):
            if nucleotide in nucleotide_map:
                one_hot[i, nucleotide_map[nucleotide]] = 1.0
        
        return one_hot
    
    def extract_sequence_with_variants(self, interval, variants):
        """Extract sequence with variants applied.
        
        Args:
            interval: Genomic interval to extract.
            variants: List of variant tuples (pos, ref, alt) to apply.
            
        Returns:
            DNA sequence as a string with variants applied.
        """
        if self.fasta_extractor is None:
            raise ValueError("Reference genome not available. Cannot extract sequence.")
        
        # Extract reference sequence
        ref_sequence = self.fasta_extractor.extract(interval)
        
        # Apply variants (shifted to be relative to the interval)
        # Sort variants by position in descending order to avoid position shifting
        variants_sorted = sorted(variants, key=lambda v: v[0], reverse=True)
        
        # Create a new sequence with variants
        variant_sequence = ref_sequence
        for var_pos, var_ref, var_alt in variants_sorted:
            # Calculate position relative to interval
            rel_pos = var_pos - interval.start
            
            # Skip variants outside the interval
            if rel_pos < 0 or rel_pos >= len(variant_sequence):
                continue
                
            # Verify reference matches
            if variant_sequence[rel_pos:rel_pos + len(var_ref)].upper() != var_ref.upper():
                print(f"Warning: Reference mismatch at position {var_pos}. "
                      f"Expected {var_ref}, found {variant_sequence[rel_pos:rel_pos + len(var_ref)]}")
                continue
            
            # Apply variant
            variant_sequence = (variant_sequence[:rel_pos] + 
                              var_alt + 
                              variant_sequence[rel_pos + len(var_ref):])
        
        return variant_sequence
    
    def predict_sequence(self, sequence):
        """Make predictions for a sequence using Enformer.
        
        Args:
            sequence: DNA sequence as a string.
            
        Returns:
            Numpy array of predictions.
        """
        # One-hot encode sequence
        one_hot = self.one_hot_encode(sequence)[np.newaxis]
        
        # Make predictions
        predictions = self.model.predict_on_batch(one_hot)['human'][0]
        
        return predictions
    
    def process_vcf_region(self, vcf_path, region, sample_limit=None):
        """Process variants from a VCF file in a specific region.
        
        Args:
            vcf_path: Path to the VCF file.
            region: Genomic region as string "chr:start-end".
            sample_limit: Maximum number of samples to process (None for all).
            
        Returns:
            Dictionary with samples and their variants.
        """
        # Parse the region
        chrom, pos_range = region.split(':')
        start, end = map(int, pos_range.split('-'))
        
        # Create a genomic interval
        interval = Interval(chrom, start, end)
        
        # Open the VCF file
        vcf_reader = vcf.Reader(filename=vcf_path)
        
        # Extract variants in the region
        variants_by_sample = {}
        
        try:
            # Try to use the built-in fetch method if the VCF is indexed
            records = vcf_reader.fetch(chrom, start, end)
        except:
            # If not indexed, filter manually (less efficient)
            print("VCF file is not indexed. Filtering manually (this may be slow).")
            records = [record for record in vcf_reader 
                      if record.CHROM == chrom and start <= record.POS <= end]
        
        # Get list of samples
        samples = vcf_reader.samples
        if sample_limit:
            samples = samples[:sample_limit]
        
        # Initialize variant dictionary for each sample
        for sample in samples:
            variants_by_sample[sample] = []
        
        # Process each variant
        for record in records:
            # For each sample, check if they have the variant
            for sample in samples:
                # Check if sample has any alternate allele
                call = record.genotype(sample)
                
                # If sample has at least one alternate allele
                if not call.is_ref:
                    # Get the alternate allele indices (could be multiple in case of heterozygous)
                    alt_indices = [i for i, gt in enumerate(call.gt_alleles) if gt != '0']
                    
                    # For simplicity, just use the first alternate allele
                    if alt_indices:
                        alt_idx = int(call.gt_alleles[alt_indices[0]]) - 1  # 0-based index
                        alt_allele = record.ALT[alt_idx]
                        
                        # Add variant to the sample's list
                        variants_by_sample[sample].append(
                            (record.POS, record.REF, str(alt_allele))
                        )
        
        return {
            'interval': interval,
            'variants_by_sample': variants_by_sample
        }
    
    def analyze_multiple_samples(self, vcf_path, region, sample_limit=None, 
                                tracks_to_analyze=None):
        """Analyze multiple samples from a VCF file.
        
        Args:
            vcf_path: Path to the VCF file.
            region: Genomic region as string "chr:start-end".
            sample_limit: Maximum number of samples to process (None for all).
            tracks_to_analyze: List of track indices to analyze. If None, a default set is used.
            
        Returns:
            Dictionary with analysis results.
        """
        if self.fasta_extractor is None:
            raise ValueError("Reference genome not available. Cannot analyze variants.")
        
        # Process VCF to get variants by sample
        vcf_data = self.process_vcf_region(vcf_path, region, sample_limit)
        interval = vcf_data['interval']
        variants_by_sample = vcf_data['variants_by_sample']
        
        # Create a larger interval for Enformer input
        # Enformer needs a larger context to make predictions
        enformer_interval = interval.resize(SEQUENCE_LENGTH)
        
        # Check if we need to select tracks
        if tracks_to_analyze is None:
            # Find some interesting tracks
            tracks_to_analyze = []
            
            # Add some CAGE tracks (gene expression)
            cage_indices = self.targets_df[self.targets_df['type'] == 'CAGE'].iloc[:5]['index']
            tracks_to_analyze.extend(cage_indices.astype(int).tolist())
            
            # Add some DNASE tracks (chromatin accessibility)
            dnase_indices = self.targets_df[self.targets_df['type'] == 'DNASE'].iloc[:5]['index']
            tracks_to_analyze.extend(dnase_indices.astype(int).tolist())
            
            # Add some histone mark tracks
            histone_marks = ['H3K27ac', 'H3K4me3', 'H3K4me1']
            for mark in histone_marks:
                mask = self.targets_df['description'].str.contains(mark, case=False)
                if mask.any():
                    indices = self.targets_df[mask].iloc[:2]['index'].astype(int).tolist()
                    tracks_to_analyze.extend(indices)
        
        print(f"Analyzing {len(tracks_to_analyze)} tracks for {len(variants_by_sample)} samples")
        
        # Get reference sequence and predict
        ref_sequence = self.fasta_extractor.extract(enformer_interval)
        ref_predictions = self.predict_sequence(ref_sequence)
        
        # Store predictions for each sample
        sample_predictions = {}
        diff_predictions = {}
        
        # Process each sample
        for sample, variants in variants_by_sample.items():
            print(f"Processing sample {sample} with {len(variants)} variants...")
            
            # Extract sequence with variants
            sample_sequence = self.extract_sequence_with_variants(enformer_interval, variants)
            
            # Make predictions
            sample_predictions[sample] = self.predict_sequence(sample_sequence)
            
            # Calculate differences from reference
            diff_predictions[sample] = sample_predictions[sample] - ref_predictions
        
        # Create a features matrix for analysis
        # For each sample, we'll take the maximum absolute difference for each track
        feature_matrix = np.zeros((len(variants_by_sample), len(tracks_to_analyze)))
        
        for i, sample in enumerate(variants_by_sample.keys()):
            for j, track_idx in enumerate(tracks_to_analyze):
                # Take the maximum absolute difference for this track
                feature_matrix[i, j] = np.abs(diff_predictions[sample][:, track_idx]).max()
        
        # Create a DataFrame for the feature matrix
        feature_df = pd.DataFrame(
            feature_matrix, 
            index=list(variants_by_sample.keys()),
            columns=[f"Track_{idx}" for idx in tracks_to_analyze]
        )
        
        # Add track descriptions to the results
        track_descriptions = {}
        for track_idx in tracks_to_analyze:
            track_descriptions[track_idx] = self.targets_df.iloc[track_idx]['description']
        
        # Perform statistical analysis to identify significant changes
        significant_tracks = self._identify_significant_tracks(feature_matrix, tracks_to_analyze)
        
        # Cluster samples based on regulatory profiles
        clusters = self._cluster_samples(feature_matrix)
        
        # Prepare results
        results = {
            'region': region,
            'interval': interval,
            'enformer_interval': enformer_interval,
            'reference_predictions': ref_predictions,
            'sample_predictions': sample_predictions,
            'diff_predictions': diff_predictions,
            'feature_matrix': feature_matrix,
            'feature_df': feature_df,
            'tracks_analyzed': tracks_to_analyze,
            'track_descriptions': track_descriptions,
            'significant_tracks': significant_tracks,
            'clusters': clusters,
            'samples': list(variants_by_sample.keys()),
            'variants_by_sample': variants_by_sample
        }
        
        return results
    
    def _identify_significant_tracks(self, feature_matrix, track_indices, threshold=0.05):
        """Identify tracks with significant changes across samples.
        
        Args:
            feature_matrix: Matrix of features (samples x tracks).
            track_indices: List of track indices.
            threshold: P-value threshold for significance.
            
        Returns:
            List of significant tracks with statistics.
        """
        significant_tracks = []
        
        # For each track, perform a one-sample t-test against 0
        for j, track_idx in enumerate(track_indices):
            # Extract feature values for this track across all samples
            track_values = feature_matrix[:, j]
            
            # Perform t-test against 0
            t_stat, p_value = stats.ttest_1samp(track_values, 0)
            
            # If significant, add to the list
            if p_value < threshold:
                significant_tracks.append({
                    'track_idx': track_idx,
                    'description': self.targets_df.iloc[track_idx]['description'],
                    't_statistic': t_stat,
                    'p_value': p_value,
                    'mean': np.mean(track_values),
                    'std': np.std(track_values)
                })
        
        # Sort by p-value
        significant_tracks.sort(key=lambda x: x['p_value'])
        
        return significant_tracks
    
    def _cluster_samples(self, feature_matrix, n_clusters=2):
        """Cluster samples based on their regulatory profiles.
        
        Args:
            feature_matrix: Matrix of features (samples x tracks).
            n_clusters: Number of clusters for K-means.
            
        Returns:
            Dictionary with clustering results.
        """
        # Normalize the feature matrix
        from sklearn.preprocessing import StandardScaler
        scaler = StandardScaler()
        normalized_features = scaler.fit_transform(feature_matrix)
        
        # Perform PCA for visualization
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(normalized_features)
        
        # Perform t-SNE for better visualization of complex relationships
        tsne = TSNE(n_components=2, random_state=42)
        tsne_result = tsne.fit_transform(normalized_features)
        
        # Perform K-means clustering
        kmeans = KMeans(n_clusters=n_clusters, random_state=42)
        cluster_labels = kmeans.fit_predict(normalized_features)
        
        # Prepare results
        results = {
            'pca': pca_result,
            'tsne': tsne_result,
            'cluster_labels': cluster_labels,
            'n_clusters': n_clusters,
            'pca_explained_variance': pca.explained_variance_ratio_
        }
        
        return results
    
    def identify_biomarkers(self, results, top_n=10):
        """Identify potential biomarkers from the analysis results.
        
        Args:
            results: Results from analyze_multiple_samples.
            top_n: Number of top biomarkers to return.
            
        Returns:
            DataFrame with biomarker information.
        """
        # Start with significant tracks
        significant_tracks = results['significant_tracks']
        
        if not significant_tracks:
            print("No significant tracks found.")
            return pd.DataFrame()
        
        # Get cluster information
        clusters = results['clusters']
        cluster_labels = clusters['cluster_labels']
        
        # For each track, calculate cluster separation
        for track in significant_tracks:
            track_idx = track['track_idx']
            track_values = results['feature_matrix'][:, results['tracks_analyzed'].index(track_idx)]
            
            # Calculate cluster means
            cluster_means = []
            for c in range(clusters['n_clusters']):
                cluster_samples = cluster_labels == c
                cluster_means.append(np.mean(track_values[cluster_samples]))
            
            # Calculate cluster separation (max difference between means)
            track['cluster_separation'] = max(cluster_means) - min(cluster_means)
        
        # Sort by a combined score of significance and cluster separation
        for track in significant_tracks:
            # Combine p-value and cluster separation into a single score
            # We take -log10(p_value) so that more significant tracks have higher scores
            track['biomarker_score'] = -np.log10(track['p_value']) * track['cluster_separation']
        
        # Sort by biomarker score
        significant_tracks.sort(key=lambda x: x['biomarker_score'], reverse=True)
        
        # Convert to DataFrame
        biomarkers_df = pd.DataFrame(significant_tracks[:top_n])
        
        return biomarkers_df
    
    def plot_results(self, results, output_dir='enformer_results'):
        """Plot the results of the analysis.
        
        Args:
            results: Results from analyze_multiple_samples.
            output_dir: Directory to save the plots.
        """
        os.makedirs(output_dir, exist_ok=True)
        
        # 1. Plot sample clustering
        self._plot_sample_clustering(results, os.path.join(output_dir, 'sample_clustering.png'))
        
        # 2. Plot top biomarkers
        biomarkers = self.identify_biomarkers(results)
        if not biomarkers.empty:
            self._plot_biomarkers(results, biomarkers, os.path.join(output_dir, 'biomarkers.png'))
        
        # 3. Plot regulatory profiles
        self._plot_regulatory_profiles(results, os.path.join(output_dir, 'regulatory_profiles.png'))
        
        # 4. Plot variant effects for a subset of samples
        sample_subset = results['samples'][:min(3, len(results['samples']))]
        for sample in sample_subset:
            self._plot_sample_effects(
                results, sample, 
                os.path.join(output_dir, f'sample_{sample}_effects.png')
            )
        
        print(f"Plots saved to {output_dir}")
    
    def _plot_sample_clustering(self, results, output_file):
        """Plot sample clustering results.
        
        Args:
            results: Results from analyze_multiple_samples.
            output_file: Path to save the plot.
        """
        clusters = results['clusters']
        samples = results['samples']
        
        # Create a figure with 2 subplots (PCA and t-SNE)
        fig, axes = plt.subplots(1, 2, figsize=(16, 7))
        
        # Plot PCA
        scatter = axes[0].scatter(
            clusters['pca'][:, 0], 
            clusters['pca'][:, 1], 
            c=clusters['cluster_labels'], 
            cmap='viridis', 
            s=80, 
            alpha=0.8
        )
        axes[0].set_title('PCA of Regulatory Profiles')
        axes[0].set_xlabel(f'PC1 ({clusters["pca_explained_variance"][0]:.2%} variance)')
        axes[0].set_ylabel(f'PC2 ({clusters["pca_explained_variance"][1]:.2%} variance)')
        
        # Add sample names to PCA plot
        for i, sample in enumerate(samples):
            axes[0].annotate(
                sample, 
                (clusters['pca'][i, 0], clusters['pca'][i, 1]),
                fontsize=8
            )
        
        # Plot t-SNE
        scatter = axes[1].scatter(
            clusters['tsne'][:, 0], 
            clusters['tsne'][:, 1], 
            c=clusters['cluster_labels'], 
            cmap='viridis', 
            s=80, 
            alpha=0.8
        )
        axes[1].set_title('t-SNE of Regulatory Profiles')
        axes[1].set_xlabel('t-SNE 1')
        axes[1].set_ylabel('t-SNE 2')
        
        # Add sample names to t-SNE plot
        for i, sample in enumerate(samples):
            axes[1].annotate(
                sample, 
                (clusters['tsne'][i, 0], clusters['tsne'][i, 1]),
                fontsize=8
            )
        
        # Add a color bar
        cbar = plt.colorbar(scatter, ax=axes.ravel().tolist())
        cbar.set_label('Cluster')
        
        plt.suptitle('Clustering of Samples Based on Regulatory Profiles', fontsize=16)
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        
        plt.savefig(output_file)
        plt.close()
    
    def _plot_biomarkers(self, results, biomarkers, output_file):
        """Plot biomarker information.
        
        Args:
            results: Results from analyze_multiple_samples.
            biomarkers: DataFrame with biomarker information.
            output_file: Path to save the plot.
        """
        # Plot biomarker scores
        plt.figure(figsize=(12, 8))
        
        # Create a bar plot of biomarker scores
        sns.barplot(
            x='biomarker_score', 
            y='description', 
            data=biomarkers.sort_values('biomarker_score')
        )
        
        plt.title('Top Regulatory Biomarkers', fontsize=16)
        plt.xlabel('Biomarker Score (-log10(p) × cluster separation)')
        plt.ylabel('Track Description')
        plt.tight_layout()
        
        plt.savefig(output_file)
        plt.close()
        
        # Create a heatmap of top biomarkers across samples
        plt.figure(figsize=(14, 10))
        
        # Get the top biomarker indices
        top_indices = [track['track_idx'] for track in biomarkers.to_dict('records')]
        
        # Create a matrix for the heatmap
        heatmap_data = np.zeros((len(results['samples']), len(top_indices)))
        
        for i, sample in enumerate(results['samples']):
            for j, track_idx in enumerate(top_indices):
                # Find the track in the tracks_analyzed list
                track_pos = results['tracks_analyzed'].index(track_idx)
                heatmap_data[i, j] = results['feature_matrix'][i, track_pos]
        
        # Create a DataFrame for the heatmap
        heatmap_df = pd.DataFrame(
            heatmap_data,
            index=results['samples'],
            columns=[biomarkers.iloc[j]['description'] for j in range(len(top_indices))]
        )
        
        # Plot the heatmap
        sns.clustermap(
            heatmap_df,
            cmap='viridis',
            figsize=(14, 10),
            row_cluster=True,
            col_cluster=True,
            xticklabels=1,
            yticklabels=1
        )
        
        plt.savefig(output_file.replace('.png', '_heatmap.png'))
        plt.close()
    
    def _plot_regulatory_profiles(self, results, output_file):
        """Plot regulatory profiles for all samples.
        
        Args:
            results: Results from analyze_multiple_samples.
            output_file: Path to save the plot.
        """
        # Get the top tracks to plot
        if results['significant_tracks']:
            top_tracks = [track['track_idx'] for track in results['significant_tracks'][:3]]
        else:
            # If no significant tracks, just take the first 3
            top_tracks = results['tracks_analyzed'][:3]
        
        # Create a figure with 3 subplots (one for each track)
        fig, axes = plt.subplots(len(top_tracks), 1, figsize=(14, 4 * len(top_tracks)), sharex=True)
        
        if len(top_tracks) == 1:
            axes = [axes]
        
        # Plot each track
        for i, track_idx in enumerate(top_tracks):
            ax = axes[i]
            track_desc = results['track_descriptions'][track_idx]
            
            # Plot reference prediction
            ax.plot(
                results['reference_predictions'][:, track_idx],
                color='black',
                alpha=0.7,
                linewidth=2,
                label='Reference'
            )
            
            # Plot predictions for each sample
            for sample in results['samples']:
                ax.plot(
                    results['sample_predictions'][sample][:, track_idx],
                    alpha=0.3,
                    linewidth=1,
                    label=sample
                )
            
            ax.set_title(f'Track: {track_desc}')
            ax.set_ylabel('Predicted Value')
            
            # Add legend to the first plot only
            if i == 0:
                ax.legend(loc='upper right')
        
        axes[-1].set_xlabel('Position (bin)')
        plt.suptitle('Regulatory Profiles Across Samples', fontsize=16)
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        
        plt.savefig(output_file)
        plt.close()
    
    def _plot_sample_effects(self, results, sample, output_file):
        """Plot effects for a specific sample.
        
        Args:
            results: Results from analyze_multiple_samples.
            sample: Sample name.
            output_file: Path to save the plot.
        """
        # Get the top tracks to plot
        if results['significant_tracks']:
            top_tracks = [track['track_idx'] for track in results['significant_tracks'][:4]]
        else:
            # If no significant tracks, just take the first 4
            top_tracks = results['tracks_analyzed'][:4]
        
        # Create a figure with subplots
        fig, axes = plt.subplots(len(top_tracks), 1, figsize=(14, 3 * len(top_tracks)), sharex=True)
        
        if len(top_tracks) == 1:
            axes = [axes]
        
        # Plot differences for each track
        for i, track_idx in enumerate(top_tracks):
            ax = axes[i]
            track_desc = results['track_descriptions'][track_idx]
            
            # Plot difference
            ax.fill_between(
                range(len(results['diff_predictions'][sample])),
                results['diff_predictions'][sample][:, track_idx],
                0,
                where=(results['diff_predictions'][sample][:, track_idx] > 0),
                color='red',
                alpha=0.5,
                label='Increase'
            )
            
            ax.fill_between(
                range(len(results['diff_predictions'][sample])),
                results['diff_predictions'][sample][:, track_idx],
                0,
                where=(results['diff_predictions'][sample][:, track_idx] < 0),
                color='blue',
                alpha=0.5,
                label='Decrease'
            )
            
            ax.set_title(f'Track: {track_desc}')
            ax.set_ylabel('Difference from Reference')
            
            # Add legend to the first plot only
            if i == 0:
                ax.legend(loc='upper right')
            
            # Add horizontal line at y=0
            ax.axhline(y=0, color='black', linestyle='-', alpha=0.3)
            
            # Mark variant positions
            variants = results['variants_by_sample'][sample]
            for var_pos, _, _ in variants:
                # Convert genomic position to bin index
                bin_idx = (var_pos - results['enformer_interval'].start) // BIN_SIZE
                if 0 <= bin_idx < len(results['diff_predictions'][sample]):
                    ax.axvline(x=bin_idx, color='green', linestyle='--', alpha=0.5)
        
        axes[-1].set_xlabel('Position (bin)')
        plt.suptitle(f'Regulatory Effects for Sample: {sample}', fontsize=16)
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        
        plt.savefig(output_file)
        plt.close()

# Example usage with mock data (since we need VCF and reference genome)
def example():
    """Example showing how to use the EnformerMultiSampleAnalyzer."""
    analyzer = EnformerMultiSampleAnalyzer(
        model_path='https://tfhub.dev/deepmind/enformer/1',
        genome_fasta_path='path_to_genome.fa'  # Replace with actual path
    )
    
    # For demonstration purposes with mock data
    print("This is a demonstration. For real analysis, provide paths to:")
    print("1. Reference genome file")
    print("2. VCF file with sample variants")
    
    # Create a mock VCF reader and data
    class MockVCFReader:
        def __init__(self):
            self.samples = ['Sample1', 'Sample2', 'Sample3', 'Sample4', 'Sample5']
            
        def fetch(self, chrom, start, end):
            # Generate some mock variants
            mock_variants = []
            import random
            
            # Generate 20 random variants in the region
            for _ in range(20):
                pos = random.randint(start, end)
                ref = random.choice(['A', 'C', 'G', 'T'])
                alt = random.choice([a for a in ['A', 'C', 'G', 'T'] if a != ref])
                
                # Create a mock variant record
                mock_variants.append(MockVariant(chrom, pos, ref, alt))
            
            return mock_variants
    
    class MockVariant:
        def __init__(self, chrom, pos, ref, alt):
            self.CHROM = chrom
            self.POS = pos
            self.REF = ref
            self.ALT = [alt]
            
        def genotype(self, sample):
            # Randomly determine if this sample has this variant
            import random
            has_variant = random.choice([True, False])
            
            return MockCall(has_variant)
    
    class MockCall:
        def __init__(self, has_variant):
            self.is_ref = not has_variant
            self.gt_alleles = ['0', '1'] if has_variant else ['0', '0']
    
    # Create a mock sequence extractor
    class MockFastaExtractor:
        def extract(self, interval):
            # Generate a random sequence of the correct length
            import random
            sequence = ''.join(random.choice(['A', 'C', 'G', 'T']) for _ in range(interval.length))
            return sequence
    
    # Set up mock objects
    analyzer.fasta_extractor = MockFastaExtractor()
    
    # Mock the VCF processing
    def mock_process_vcf_region(self, vcf_path, region, sample_limit=None):
        # Parse the region
        chrom, pos_range = region.split(':')
        start, end = map(int, pos_range.split('-'))
        
        # Create a genomic interval
        interval = Interval(chrom, start, end)
        
        # Create mock samples and variants
        vcf_reader = MockVCFReader()
        samples = vcf_reader.samples[:sample_limit] if sample_limit else vcf_reader.samples
        
        # Generate random variants for each sample
        variants_by_sample = {}
        import random
        
        for sample in samples:
            variants_by_sample[sample] = []
            # Generate 3-10 random variants for this sample
            for _ in range(random.randint(3, 10)):
                pos = random.randint(start, end)
                ref = random.choice(['A', 'C', 'G', 'T'])
                alt = random.choice([a for a in ['A', 'C', 'G', 'T'] if a != ref])
                variants_by_sample[sample].append((pos, ref, alt))
        
        return {
            'interval': interval,
            'variants_by_sample': variants_by_sample
        }
    
    # Replace the real method with our mock version for the demo
    analyzer.process_vcf_region = mock_process_vcf_region.__get__(analyzer)
    
    # Also mock the predict_sequence method to avoid actual model inference
    def mock_predict_sequence(self, sequence):
        # Generate random predictions
        import numpy as np
        # Generate predictions of shape [896, 5313]
        predictions = np.random.random((896, 5313)) * 0.1
        return predictions
    
    # Replace the real method with our mock version for the demo
    analyzer.predict_sequence = mock_predict_sequence.__get__(analyzer)
    
    # Run the analysis on a mock region
    results = analyzer.analyze_multiple_samples(
        vcf_path='mock.vcf',
        region='chr7:55,000,000-55,050,000',  # EGFR gene region (mock)
        sample_limit=5
    )
    
    # Identify potential biomarkers
    biomarkers = analyzer.identify_biomarkers(results)
    print("Top biomarkers:")
    if not biomarkers.empty:
        print(biomarkers[['description', 'biomarker_score', 'p_value']].head())
    
    # Plot the results
    analyzer.plot_results(results, 'enformer_results')
    
    print("\nAnalysis complete!")
    print("In a real-world scenario, this pipeline would:")
    print("1. Process real VCF files with patient/sample variants")
    print("2. Use Enformer to predict regulatory effects for each sample")
    print("3. Identify biomarkers that distinguish between sample groups")
    print("4. Cluster samples based on their regulatory profiles")
    print("5. Generate visualizations to interpret the results")

if __name__ == "__main__":
    example()




