"""
Enformer Variant Analysis: Detailed Example

This script demonstrates how to analyze genetic variants using Enformer by:
1. Creating reference and variant sequences
2. Running both through Enformer
3. Comparing predictions to quantify variant effects
4. Visualizing the differences
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tensorflow as tf
import tensorflow_hub as hub
from kipoiseq import Interval
from kipoiseq.extractors import FastaStringExtractor
from kipoiseq.transforms.functional import one_hot_encode
import seaborn as sns
from pyfaidx import Fasta

# Constants
SEQUENCE_LENGTH = 393_216  # Input sequence length for Enformer
BIN_SIZE = 128  # Size of bins for prediction

class EnformerVariantAnalyzer:
    """Class for analyzing genetic variants using Enformer."""
    
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
        
        # Check if file exists before loading
        if os.path.exists(genome_fasta_path):
            print(f"Loading reference genome from {genome_fasta_path}...")
            self.fasta_extractor = FastaStringExtractor(genome_fasta_path)
        else:
            print(f"Warning: Genome file {genome_fasta_path} not found.")
            self.fasta_extractor = None
            
        # Load target information if available
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
    
    def _create_variant_sequence(self, sequence, variant_position, ref_allele, alt_allele):
        """Create a sequence with the variant.
        
        Args:
            sequence: Reference DNA sequence.
            variant_position: Position in the sequence to introduce the variant (0-based).
            ref_allele: Reference allele.
            alt_allele: Alternate allele.
            
        Returns:
            Sequence with the variant.
        """
        # Verify the reference allele matches what's in the sequence
        seq_ref = sequence[variant_position:variant_position + len(ref_allele)]
        if seq_ref.upper() != ref_allele.upper():
            raise ValueError(f"Reference allele mismatch at position {variant_position}. "
                             f"Expected {ref_allele}, found {seq_ref}.")
        
        # Create variant sequence
        return sequence[:variant_position] + alt_allele + sequence[variant_position + len(ref_allele):]
    
    def analyze_variant(self, chrom, pos, ref_allele, alt_allele, 
                        window_size=200_000, tracks_to_analyze=None):
        """Analyze a variant using Enformer.
        
        Args:
            chrom: Chromosome name (e.g., 'chr1').
            pos: Genomic position (1-based).
            ref_allele: Reference allele.
            alt_allele: Alternate allele.
            window_size: Size of the window to extract around the variant.
            tracks_to_analyze: List of track indices to analyze. If None, a default set is used.
            
        Returns:
            Dictionary with analysis results.
        """
        if self.fasta_extractor is None:
            raise ValueError("Reference genome not available. Cannot analyze variant.")
        
        # Create interval centered on the variant
        variant_interval = Interval(chrom, pos - window_size//2, pos + window_size//2)
        variant_centered = variant_interval.resize(SEQUENCE_LENGTH)
        
        # Extract reference sequence
        ref_sequence = self.fasta_extractor.extract(variant_centered)
        
        # Calculate position in the extracted sequence
        variant_position = pos - variant_centered.start
        
        # Create variant sequence
        alt_sequence = self._create_variant_sequence(
            ref_sequence, variant_position, ref_allele, alt_allele)
        
        # One-hot encode sequences
        ref_one_hot = one_hot_encode(ref_sequence)[np.newaxis]
        alt_one_hot = one_hot_encode(alt_sequence)[np.newaxis]
        
        # Make predictions
        print("Making predictions for reference sequence...")
        ref_predictions = self.model.predict_on_batch(ref_one_hot)['human'][0]
        
        print("Making predictions for alternate sequence...")
        alt_predictions = self.model.predict_on_batch(alt_one_hot)['human'][0]
        
        # Calculate differences
        diff_predictions = alt_predictions - ref_predictions
        
        # Select tracks to analyze if not specified
        if tracks_to_analyze is None:
            # Find some interesting tracks
            tracks_to_analyze = []
            
            # Add a CAGE track (gene expression)
            cage_idx = self.targets_df[self.targets_df['type'] == 'CAGE'].iloc[0]['index']
            tracks_to_analyze.append(int(cage_idx))
            
            # Add a DNASE track (chromatin accessibility)
            dnase_idx = self.targets_df[self.targets_df['type'] == 'DNASE'].iloc[0]['index']
            tracks_to_analyze.append(int(dnase_idx))
            
            # Add some histone mark tracks if available
            histone_marks = ['H3K27ac', 'H3K4me3', 'H3K4me1']
            for mark in histone_marks:
                mask = self.targets_df['description'].str.contains(mark, case=False)
                if mask.any():
                    tracks_to_analyze.append(int(self.targets_df[mask].iloc[0]['index']))
        
        # Prepare results
        results = {
            'variant': f"{chrom}:{pos} {ref_allele}>{alt_allele}",
            'ref_predictions': ref_predictions,
            'alt_predictions': alt_predictions,
            'diff_predictions': diff_predictions,
            'tracks_analyzed': tracks_to_analyze,
            'track_info': {idx: self.targets_df.iloc[idx]['description'] for idx in tracks_to_analyze},
            'variant_position': variant_position,
            'variant_interval': variant_centered,
        }
        
        # Calculate summary statistics
        results['max_abs_diff'] = np.abs(diff_predictions).max(axis=0)
        results['mean_abs_diff'] = np.abs(diff_predictions).mean(axis=0)
        
        # Calculate effect size for the analyzed tracks
        results['track_effects'] = {}
        for track_idx in tracks_to_analyze:
            # Find the bin with the largest change for this track
            max_change_bin = np.argmax(np.abs(diff_predictions[:, track_idx]))
            max_change_value = diff_predictions[max_change_bin, track_idx]
            
            # Calculate relative change
            baseline = ref_predictions[max_change_bin, track_idx]
            if baseline != 0:
                relative_change = max_change_value / baseline
            else:
                relative_change = float('inf') if max_change_value > 0 else float('-inf')
            
            results['track_effects'][track_idx] = {
                'max_change_bin': max_change_bin,
                'max_change_value': max_change_value,
                'relative_change': relative_change,
                'bin_position': variant_centered.start + max_change_bin * BIN_SIZE
            }
        
        return results
    
    def plot_variant_effects(self, results, output_file=None):
        """Plot the effects of a variant.
        
        Args:
            results: Results from analyze_variant.
            output_file: Path to save the plot. If None, plot is displayed.
        """
        tracks_to_plot = results['tracks_analyzed']
        n_tracks = len(tracks_to_plot)
        
        # Create a figure with three columns (ref, alt, diff) and n_tracks rows
        fig, axes = plt.subplots(n_tracks, 3, figsize=(18, 4 * n_tracks))
        
        # If only one track, make axes 2D
        if n_tracks == 1:
            axes = axes.reshape(1, -1)
        
        # Plot each track
        for i, track_idx in enumerate(tracks_to_plot):
            track_name = results['track_info'][track_idx]
            
            # Plot reference
            axes[i, 0].fill_between(
                range(results['ref_predictions'].shape[0]),
                results['ref_predictions'][:, track_idx]
            )
            axes[i, 0].set_title(f"Reference: {track_name}")
            
            # Plot alternate
            axes[i, 1].fill_between(
                range(results['alt_predictions'].shape[0]),
                results['alt_predictions'][:, track_idx]
            )
            axes[i, 1].set_title(f"Alternate: {track_name}")
            
            # Plot difference
            axes[i, 2].fill_between(
                range(results['diff_predictions'].shape[0]),
                results['diff_predictions'][:, track_idx]
            )
            axes[i, 2].set_title(f"Difference: {track_name}")
            
            # Mark the variant position
            variant_bin = results['variant_position'] // BIN_SIZE
            for j in range(3):
                axes[i, j].axvline(x=variant_bin, color='red', linestyle='--')
                
            # Annotate the effect size
            effect = results['track_effects'][track_idx]
            axes[i, 2].annotate(
                f"Max change: {effect['max_change_value']:.4f}\n"
                f"Relative change: {effect['relative_change']:.2%}",
                xy=(0.05, 0.95),
                xycoords='axes fraction',
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
            )
        
        plt.suptitle(f"Variant Effects Analysis: {results['variant']}", fontsize=16)
        plt.tight_layout(rect=[0, 0, 1, 0.97])
        
        if output_file:
            plt.savefig(output_file)
            print(f"Plot saved to {output_file}")
        else:
            plt.show()
    
    def analyze_multiple_variants(self, vcf_path, region=None, max_variants=10):
        """Analyze multiple variants from a VCF file.
        
        Args:
            vcf_path: Path to the VCF file.
            region: Genomic region to analyze (e.g., 'chr1:1000000-2000000').
            max_variants: Maximum number of variants to analyze.
            
        Returns:
            DataFrame with analysis results.
        """
        # This is a placeholder for VCF processing functionality
        # In a real implementation, you would use a library like PyVCF
        # to parse the VCF file and analyze each variant
        print(f"Analyzing variants from {vcf_path}")
        print("This functionality requires implementing a VCF parser.")
        print("For a full implementation, consider using the PyVCF library.")
        
        # Return mock results
        return pd.DataFrame({
            'variant': [f"chr1:{1000000+i} A>G" for i in range(max_variants)],
            'max_effect': np.random.random(max_variants),
            'affected_track': ["CAGE_Brain"] * max_variants
        })

# Example usage
def example():
    """Example of using the EnformerVariantAnalyzer."""
    analyzer = EnformerVariantAnalyzer(
        model_path='https://tfhub.dev/deepmind/enformer/1',
        genome_fasta_path='hg38.fa'  # Replace with actual path
    )
    
    # Example variant: rs1234 (fictional)
    chrom = 'chr7'
    pos = 55_200_000  # Near EGFR gene
    ref_allele = 'A'
    alt_allele = 'G'
    
    # Mock sequence for demonstration
    # In a real scenario, you would use the actual genome sequence
    mock_sequence = ''.join(np.random.choice(['A', 'C', 'G', 'T'], size=SEQUENCE_LENGTH))
    mock_sequence = mock_sequence[:pos-55_100_000] + ref_allele + mock_sequence[pos-55_100_000+1:]
    
    # Mock implementation for demonstration
    class MockExtractor:
        def extract(self, interval):
            return mock_sequence
    
    analyzer.fasta_extractor = MockExtractor()
    
    # Analyze the variant
    results = analyzer.analyze_variant(chrom, pos, ref_allele, alt_allele)
    
    # Plot the results
    analyzer.plot_variant_effects(results, output_file='variant_effects.png')
    
    # Print summary of effects
    print(f"\nVariant {results['variant']} analysis complete:")
    for track_idx, effect in results['track_effects'].items():
        track_name = results['track_info'][track_idx]
        print(f"  {track_name}:")
        print(f"    Max change: {effect['max_change_value']:.4f}")
        print(f"    Relative change: {effect['relative_change']:.2%}")
        print(f"    Position of max effect: {effect['bin_position']}")
    
    print("\nNote: This is a demonstration with mock data.")
    print("For real analysis, provide the path to an actual reference genome.")

if __name__ == "__main__":
    example()


