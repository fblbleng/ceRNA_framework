# Optionally, filter interactions using a global threshold.
global_threshold = merged_interactions['confidence_score'].quantile(0.10)
final_df = merged_interactions[merged_interactions['confidence_score'] >= global_threshold]

# Save the final merged interaction data.
final_df.to_csv('data/final_ceRNA_interactions.csv', index=False)
