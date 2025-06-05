import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import confusion_matrix, classification_report

def evaluate_and_compare(knn_preds, wknn_preds, label_paths, output_dir='./output'):
    """评估并比较KNN和WKNN的性能"""
    # 加载真实标签
    if not label_paths:
        print("没有提供标签，无法进行评估")
        return None
    
    atac_labels = np.loadtxt(label_paths[0])
    for i in range(1, len(label_paths)):
        atac_labels = np.concatenate((atac_labels, np.loadtxt(label_paths[i])), 0)
    
    valid_indices = np.where(atac_labels >= 0)[0]
    if len(valid_indices) == 0:
        print("没有有效标签用于评估")
        return None
    
    # 计算准确率
    knn_acc = np.mean(knn_preds[valid_indices] == atac_labels[valid_indices])
    wknn_acc = np.mean(wknn_preds[valid_indices] == atac_labels[valid_indices])
    
    # 生成详细分类报告
    knn_report = classification_report(atac_labels[valid_indices], knn_preds[valid_indices])
    wknn_report = classification_report(atac_labels[valid_indices], wknn_preds[valid_indices])
    
    # 保存报告
    with open(f'{output_dir}/knn_classification_report.txt', 'w') as f:
        f.write(knn_report)
    with open(f'{output_dir}/wknn_classification_report.txt', 'w') as f:
        f.write(wknn_report)
    
    # 可视化混淆矩阵
    plt.figure(figsize=(16, 6))
    
    plt.subplot(1, 2, 1)
    cm_knn = confusion_matrix(atac_labels[valid_indices], knn_preds[valid_indices])
    sns.heatmap(cm_knn, annot=True, fmt='d', cmap='Blues')
    plt.title('KNN Confusion Matrix')
    plt.xlabel('Predicted')
    plt.ylabel('True')
    
    plt.subplot(1, 2, 2)
    cm_wknn = confusion_matrix(atac_labels[valid_indices], wknn_preds[valid_indices])
    sns.heatmap(cm_wknn, annot=True, fmt='d', cmap='Greens')
    plt.title('WKNN Confusion Matrix')
    plt.xlabel('Predicted')
    plt.ylabel('True')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/confusion_matrices_comparison.png')
    
    # 性能提升分析
    improvement = (wknn_acc - knn_acc) / knn_acc * 100
    
    return {
        'knn_accuracy': knn_acc,
        'wknn_accuracy': wknn_acc,
        'improvement_percentage': improvement
    }

def visualize_tuning_results(results_file='./output/wknn_tuning_results.csv'):
    """可视化参数调优结果"""
    if not os.path.exists(results_file):
        print(f"找不到调优结果文件: {results_file}")
        return
    
    results_df = pd.read_csv(results_file)
    
    # 按weight_type分组绘制热力图
    for wt in results_df['weight_type'].unique():
        wt_df = results_df[results_df['weight_type'] == wt]
        
        if wt == 'gaussian':
            pivot = wt_df.pivot(index='sigma', columns='k', values='accuracy')
            plt.figure(figsize=(10, 8))
            sns.heatmap(pivot, annot=True, cmap='YlGnBu', fmt='.4f')
            plt.title(f'Gaussian Weight - Accuracy by K and Sigma')
            plt.savefig(f'./output/wknn_tuning_gaussian.png')
        else:
            pivot = wt_df.pivot(index='k', columns='weight_type', values='accuracy')
            plt.figure(figsize=(10, 8))
            sns.heatmap(pivot, annot=True, cmap='YlGnBu', fmt='.4f')
            plt.title(f'Distance Weight - Accuracy by K')
            plt.savefig(f'./output/wknn_tuning_distance.png')
    
    print("参数调优可视化已保存至 ./output/")