from sklearn.neighbors import NearestNeighbors
import numpy as np
from scipy.linalg import norm
from scipy.special import softmax
import sys
import os

from config import Config
from sklearn.preprocessing import normalize

def neighbor_hit_cnt(rna_cnt, neighbor_indexs):
    hit_cnt = np.zeros(rna_cnt)
    for i in range(neighbor_indexs.shape[0]):
        for j in range(neighbor_indexs.shape[1]):
            hit_cnt[neighbor_indexs[i][j]] += 1

    return hit_cnt

def weighted_knn_predict(atac_embeddings, rna_embedding_knn, rna_label_knn, k=40):
    neigh = NearestNeighbors(n_neighbors=k)
    neigh.fit(rna_embedding_knn)
    distances, indices = neigh.kneighbors(atac_embeddings)

    predictions = []
    for i in range(len(atac_embeddings)):
        class_weights = {}
        for j, idx in enumerate(indices[i]):
            label = rna_label_knn[idx]
            dist = distances[i][j]
            weight = 1 / (dist + 1e-5)
            class_weights[label] = class_weights.get(label, 0) + weight
        predictions.append(max(class_weights, key=class_weights.get))

    return np.array(predictions), indices

def compute_hit_conf(knn_label, rna_labels, neighbor_indexs, hit_cnts):
    num_samples = knn_label.shape[0]
    conf_scores = np.zeros(num_samples)
    for i in range(num_samples):
        for j in range(neighbor_indexs.shape[1]):
            if knn_label[i] == rna_labels[neighbor_indexs[i][j]]:
                conf_scores[i] += 1/hit_cnts[neighbor_indexs[i][j]]
            else:
                conf_scores[i] -= 1/hit_cnts[neighbor_indexs[i][j]]
    return conf_scores

def KNN(config, neighbors = 30, knn_rna_samples = 20000):    
    print('[WKNN] Read RNA data')
    db_name = os.path.basename(config.rna_paths[0]).split('.')[0]
    rna_embeddings = np.loadtxt('./output/' + db_name + '_embeddings.txt')
    rna_predictions = np.loadtxt('./output/' + db_name + '_predictions.txt')
    rna_labels = np.loadtxt(config.rna_labels[0])    
    for i in range(1, len(config.rna_paths)):
        db_name = os.path.basename(config.rna_paths[i]).split('.')[0]
        rna_embeddings = np.concatenate((rna_embeddings, np.loadtxt('./output/' + db_name + '_embeddings.txt')), 0)
        rna_predictions = np.concatenate((rna_predictions, np.loadtxt('./output/' + db_name + '_predictions.txt')), 0)
        rna_labels = np.concatenate((rna_labels, np.loadtxt(config.rna_labels[i])), 0)

    rna_embedding_knn = []
    rna_label_knn = []
    rna_prediction_knn = []
    
    num_of_rna = rna_embeddings.shape[0]
    if num_of_rna > knn_rna_samples:
        sampling_interval = num_of_rna*1./knn_rna_samples
        i = 0
        while i < num_of_rna:        
            rna_embedding_knn.append(rna_embeddings[i])
            rna_label_knn.append(int(rna_labels[i]))
            rna_prediction_knn.append(rna_predictions[i])
            i = int(i + sampling_interval)
    else: 
        rna_embedding_knn = rna_embeddings
        rna_label_knn = rna_labels
        rna_prediction_knn = rna_predictions

    print('[WKNN] Read ATAC data')
    db_names = []
    db_sizes = []
    db_name = os.path.basename(config.atac_paths[0]).split('.')[0]    
    atac_embeddings = np.loadtxt('./output/' + db_name + '_embeddings.txt')
    atac_predictions = np.loadtxt('./output/' + db_name + '_predictions.txt')
    db_names.append(db_name)
    db_sizes.append(atac_embeddings.shape[0])
    for i in range(1, len(config.atac_paths)):
        db_name = os.path.basename(config.atac_paths[i]).split('.')[0]        
        em = np.loadtxt('./output/' + db_name + '_embeddings.txt')
        pred = np.loadtxt('./output/' + db_name + '_predictions.txt')
        atac_embeddings = np.concatenate((atac_embeddings, em), 0)
        atac_predictions = np.concatenate((atac_predictions, pred), 0)        
        db_names.append(db_name)
        db_sizes.append(em.shape[0])
    
    rna_embedding_knn = normalize(np.array(rna_embedding_knn), axis=1)
    atac_embeddings = normalize(atac_embeddings, axis=1)


    print('[WKNN] Running weighted KNN')
    atac_predict, top_neighbors = weighted_knn_predict(atac_embeddings, np.array(rna_embedding_knn), np.array(rna_label_knn), neighbors)
    hit_cnts = neighbor_hit_cnt(len(rna_label_knn), top_neighbors)
    conf_scores = compute_hit_conf(atac_predict, rna_label_knn, top_neighbors, hit_cnts)

    cnt = 0
    for i, db_name in enumerate(db_names):        
        np.savetxt('./output/' + db_name + '_knn_predictions.txt', atac_predict[cnt:cnt+db_sizes[i]])
        np.savetxt('./output/' + db_name + '_knn_probs.txt', conf_scores[cnt:cnt+db_sizes[i]])
        cnt += db_sizes[i]
    

    if len(config.atac_labels[0] and config.atac_labels) == len(config.atac_paths):
        atac_labels = np.loadtxt(config.atac_labels[0])    
        for i in range(1, len(config.atac_labels)):
            atac_labels = np.concatenate((atac_labels, np.loadtxt(config.atac_labels[i])), 0)

        valid_sample_cnt = 0
        correct = 0
        for i in range(atac_predict.shape[0]):
            if atac_labels[i] >= 0:
                valid_sample_cnt += 1
                if atac_labels[i] == atac_predict[i]:
                    correct += 1

        print('wknn accuracy:', correct*1./valid_sample_cnt)

if __name__ == "__main__":
    config = Config()
    KNN(config)
