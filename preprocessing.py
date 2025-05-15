import pandas as pd

def load_and_filter(filepath, score_threshold=None):
    # Assumes a tab-delimited or CSV format; adjust sep if needed.
    df = pd.read_csv(filepath, sep='\t')
    
    # Check for a confidence score column
    if score_threshold and 'confidence_score' in df.columns:
        df = df[df['confidence_score'] >= score_threshold]
    
    # Remove duplicates and rows with missing interaction details
    df.dropna(subset=['source', 'target'], inplace=True)
    df = df.drop_duplicates(subset=['source', 'target'])
    
    return df

# Example loading for miRTarBase data
miRTarBase_df = load_and_filter('data/miRTarBase_interactions.tsv', score_threshold=0.5)
