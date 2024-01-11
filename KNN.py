import argparse
import pandas as pd
import numpy as np

def load_data(file_path):
    # Function to load data from a given file path using pandas, the data prvided is tab delimited
    return pd.read_csv(file_path, delimiter='\t')

def euclidean_distance(point1, point2):
    # build Functions to calculate the Euclidean distance between two points (vectors)
    return np.sqrt(np.sum((point1 - point2) ** 2))

def knn(train_data, test_data, k):
    # Function to perform K-Nearest Neighbors algorithm
    predictions = []

    # Iterate over each row in the test data
    for test_index, test_row in test_data.iterrows():
        distances = [] #distances will be list of tuples, where each tuple consists of type of cancer from the training data and the Euclidean distance between training data point and the actual data point from test.

        # within - Iterate over each row in the training data
        for train_index, train_row in train_data.iterrows():
            # Calculate the Euclidean distance between the test and training data points
            dist = euclidean_distance(test_row[:-1], train_row[:-1])
            distances.append((train_row['tissue'], dist))

        # Sort distances to find the k-nearest neighbors
        distances.sort(key=lambda x: x[1])
        k_neighbors = distances[:k]
        #k_neighbors intiated as a list of tuples representing the k-nearest neighbors; in other words- each tuple contains cancer type of a training data point and the Euclidean distance from the test data point.

        # Get labels of k-nearest neighbors
        neighbor_labels = [neighbor[0] for neighbor in k_neighbors]

        # Predict the label based on the majority among the neighbors
        predicted_label = max(set(neighbor_labels), key=neighbor_labels.count)
        predictions.append(predicted_label)

    return predictions

def main():
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(description='KNN for cancer prediction')
    parser.add_argument('-r', '--train_data', required=True, help='Path to training data')
    parser.add_argument('-t', '--test_data', required=True, help='Path to test data')
    parser.add_argument('-k', '--k_value', type=int, required=True, help='Number of neighbors (k) to consider')
    parser.add_argument('-o', '--output_file', required=True, help='Path to the output file for predictions')
    args = parser.parse_args()

    # Load training and test data
    train_data = load_data(args.train_data)
    test_data = load_data(args.test_data)

    # Make predictions using KNN
    predictions = knn(train_data, test_data, args.k_value)
    
    # output -> file
    # Save predictions to the output file
    test_data['predicted_tissue'] = predictions
    test_data.to_csv(args.output_file, sep='\t', index=False)

if __name__ == '__main__':
    main()
