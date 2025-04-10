"""
Enformer Usage Example: Predicting gene expression from DNA sequence

This script demonstrates how to:
1. Set up the Enformer model from TensorFlow Hub
2. Extract a DNA sequence from a genome
3. Make predictions with Enformer
4. Visualize the predictions

Requirements:
- tensorflow (2.4.1+)
- tensorflow-hub
- kipoiseq (0.5.2+)
- pyfaidx
- numpy
- pandas
- matplotlib
- seaborn
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import tensorflow as tf
import tensorflow_hub as hub
from kipoiseq import Interval
from kipoiseq.extractors import FastaStringExtractor

# Check if GPU is available
print("GPU Available: ", tf.config.list_physical_devices('GPU'))

# Constants
SEQUENCE_LENGTH = 393_216  # Input sequence length for Enformer
BIN_SIZE = 128  # Size of bins for prediction

def one_hot_encode(sequence):
    """Convert DNA sequence to one-hot encoding.
    
    Args:
        sequence: DNA sequence as a string.
        
    Returns:
        One-hot encoded sequence as numpy array.
    """
    sequence = sequence.upper()
    # Define the mapping of nucleotides to indices: A -> 0, C -> 1, G -> 2, T -> 3
    nucleotide_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    # Initialize the one-hot encoded array
    one_hot = np.zeros((len(sequence), 4), dtype=np.float32)
    
    # Fill in the one-hot encoded array
    for i, nucleotide in enumerate(sequence):
        if nucleotide in nucleotide_map:
            one_hot[i, nucleotide_map[nucleotide]] = 1.0
    
    return one_hot

def plot_tracks(tracks, interval, height=1.5):
    """Plot prediction tracks.
    
    Args:
        tracks: Dictionary of track names and values.
        interval: Genomic interval for plotting.
        height: Height per track.
    """
    fig, axes = plt.subplots(len(tracks), 1, figsize=(20, height * len(tracks)), sharex=True)
    if len(tracks) == 1:
        axes = [axes]
        
    for ax, (title, y) in zip(axes, tracks.items()):
        ax.fill_between(np.linspace(interval.start, interval.end, num=len(y)), y)
        ax.set_title(title)
        sns.despine(top=True, right=True, bottom=True)
    
    axes[-1].set_xlabel(f"{interval.chrom}:{interval.start}-{interval.end}")
    plt.tight_layout()
    return fig

def main():
    # 1. Load the Enformer model from TensorFlow Hub
    print("Loading Enformer model...")
    enformer = hub.load('https://tfhub.dev/deepmind/enformer/1')
    print("Model loaded successfully!")

    # 2. Load genome sequence (example with hg38)
    # You'll need to download the genome reference file first
    # For this example, we assume you have the hg38.fa file already
    genome_file = 'path/to/your/genome/hg38.fa'  # Replace with your actual path
    
    # If you don't have a genome file yet, you can download a sample region
    # or use a smaller test genome for this example
    if not os.path.exists(genome_file):
        print(f"Warning: Genome file {genome_file} not found. Using mock sequence for demonstration.")
        # Create a mock sequence for demonstration
        mock_sequence = ''.join(np.random.choice(['A', 'C', 'G', 'T'], size=SEQUENCE_LENGTH))
        sequence_one_hot = one_hot_encode(mock_sequence)
    else:
        # Extract sequence from reference genome
        fasta_extractor = FastaStringExtractor(genome_file)
        
        # Define a genomic interval around the SOX9 gene (example)
        # SOX9 is a gene involved in determining sex during development
        target_interval = Interval('chr17', 70_117_000, 70_317_000)  # SOX9 gene region
        
        # Resize the interval to match Enformer's input size while keeping the gene centered
        sequence_length_centered = target_interval.resize(SEQUENCE_LENGTH)
        
        # Extract the DNA sequence
        sequence = fasta_extractor.extract(sequence_length_centered)
        sequence_one_hot = one_hot_encode(sequence)
    
    # 3. Make predictions with Enformer
    print("Making predictions...")
    # Add batch dimension
    sequence_one_hot_batch = sequence_one_hot[np.newaxis]
    
    # Get predictions for human tracks
    predictions = enformer.predict_on_batch(sequence_one_hot_batch)
    human_predictions = predictions['human'][0]  # Remove batch dimension
    
    print(f"Prediction shape: {human_predictions.shape}")
    # Shape is [896, 5313] - 896 bins of 128bp each, and 5313 human output tracks
    
    # 4. Load track metadata to know what each prediction represents
    # This can be downloaded from the Basenji2 dataset
    try:
        # Try to load targets info (you need to download this file from Basenji repository)
        targets_txt = 'https://raw.githubusercontent.com/calico/basenji/0.5/manuscripts/cross2020/targets_human.txt'
        df_targets = pd.read_csv(targets_txt, sep='\t')
        print(f"Loaded {len(df_targets)} target tracks")
    except:
        print("Could not load targets info. Creating mock targets.")
        # Create mock targets data for demonstration
        df_targets = pd.DataFrame({
            'index': range(5313),
            'description': [f'Track {i}' for i in range(5313)],
            'file': ['mock.file'] * 5313,
            'type': ['CHIP'] * 2000 + ['DNASE'] * 2000 + ['CAGE'] * 1313,
            'cell_type': ['cell_type_mock'] * 5313
        })
    
    # 5. Visualize predictions for specific tracks
    # Find example tracks from different assay types
    example_tracks = {}
    
    # Find a DNASE track
    dnase_idx = df_targets[df_targets['description'].str.contains('DNASE', case=False)].iloc[0]['index']
    example_tracks[f'DNASE: {df_targets.iloc[dnase_idx]["description"]}'] = human_predictions[:, dnase_idx]
    
    # Find a CHIP-seq H3K27ac track (marker of active enhancers)
    h3k27ac_idx = df_targets[df_targets['description'].str.contains('H3K27ac', case=False)].iloc[0]['index']
    example_tracks[f'CHIP H3K27ac: {df_targets.iloc[h3k27ac_idx]["description"]}'] = human_predictions[:, h3k27ac_idx]
    
    # Find a CAGE track (marker of transcription start sites)
    cage_idx = df_targets[df_targets['description'].str.contains('CAGE', case=False)].iloc[0]['index']
    # CAGE values are typically log-transformed 
    example_tracks[f'CAGE: {df_targets.iloc[cage_idx]["description"]}'] = np.log10(1 + human_predictions[:, cage_idx])
    
    # Create the target interval for plotting (central portion of the input)
    central_start = SEQUENCE_LENGTH // 2 - (BIN_SIZE * 896) // 2
    central_end = central_start + (BIN_SIZE * 896)
    plot_interval = Interval('chr17', 70_117_000 + central_start, 70_117_000 + central_end)
    
    # Plot the tracks
    fig = plot_tracks(example_tracks, plot_interval)
    plt.savefig('enformer_predictions.png')
    print("Predictions plotted and saved to enformer_predictions.png")
    
    # 6. Optional: Compute contribution scores to see which parts of the DNA influence the prediction
    # For example, to find regions contributing to a specific CAGE track
    print("Computing contribution scores (this can take a few minutes)...")
    target_mask = np.zeros_like(human_predictions)
    
    # Focus on the central bins for the CAGE track
    central_bins = [447, 448, 449]  # Central bins
    for idx in central_bins:
        target_mask[idx, cage_idx] = 1
    
    # Compute contribution scores using gradient x input
    # Note: This is a simplified version. For better results, use the actual implementation from Enformer
    with tf.GradientTape() as tape:
        inputs = tf.constant(sequence_one_hot_batch, dtype=tf.float32)
        tape.watch(inputs)
        predictions = enformer.predict_on_batch(inputs)
        human_preds = predictions['human']
        target_pred = tf.reduce_sum(human_preds * tf.constant(target_mask, dtype=tf.float32))
    
    input_grad = tape.gradient(target_pred, inputs).numpy()
    
    # Compute contribution scores (gradient * input)
    contribution_scores = input_grad[0] * sequence_one_hot
    
    # Sum across ACGT dimension
    contribution_scores_sum = contribution_scores.sum(axis=1)
    
    # Plot contribution scores
    plt.figure(figsize=(20, 3))
    plt.plot(contribution_scores_sum)
    plt.title(f'Contribution scores for {df_targets.iloc[cage_idx]["description"]}')
    plt.xlabel('Position in sequence')
    plt.ylabel('Contribution score')
    plt.savefig('enformer_contribution_scores.png')
    print("Contribution scores plotted and saved to enformer_contribution_scores.png")
    
    print("Done!")

if __name__ == "__main__":
    main()



