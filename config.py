import torch
import os

class Config(object):
    def __init__(self):
        DB = 'HematopoiesisData'
        self.use_cuda = False
        self.threads = 1
        self.temperature = 0.1 
        self.weight_decay = 1e-4
        self.num_res_blocks = 2  # 每个编码器/解码器中的残差块数量
        self.dropout = 0.2       # Dropout率
        self.hidden_size = 512  
        if not self.use_cuda:
            self.device = torch.device('cpu')
        else:
            self.device = torch.device('cuda:0')
        
        if DB == '10x':
            # DB info
            self.number_of_class = 11
            self.input_size = 15463
            self.rna_paths = ['data_10x/exprs_10xPBMC_rna.npz']
            self.rna_labels = ['data_10x/cellType_10xPBMC_rna.txt']		
            self.atac_paths = ['data_10x/exprs_10xPBMC_atac.npz']
            self.atac_labels = [''] #Optional. If atac_labels are provided, accuracy after knn would be provided.
            self.rna_protein_paths = []
            self.atac_protein_paths = []
            
            # Training config            
            self.batch_size = 256
            self.lr_stage1 = 0.01
            self.lr_stage3 = 0.01
            self.lr_decay_epoch = 20
            self.epochs_stage1 = 20
            self.epochs_stage3 = 20
            self.p = 0.8
            self.embedding_size = 64
            self.momentum = 0.9
            self.center_weight = 1
            self.with_crossentorpy = True
            self.seed = 1
            self.checkpoint = ''
        
        elif DB == "MOp":
            self.number_of_class = 21
            self.input_size = 18603
            self.rna_paths = ['data_MOp/YaoEtAl_RNA_snRNA_10X_v3_B_exprs.npz',\
                                'data_MOp/YaoEtAl_RNA_snRNA_10X_v3_A_exprs.npz',\
                                'data_MOp/YaoEtAl_RNA_snRNA_10X_v2_exprs.npz',\
                                'data_MOp/YaoEtAl_RNA_snRNA_SMARTer_exprs.npz',\
                                'data_MOp/YaoEtAl_RNA_scRNA_10X_v3_exprs.npz',\
                                'data_MOp/YaoEtAl_RNA_scRNA_10X_v2_exprs.npz',\
                                'data_MOp/YaoEtAl_RNA_scRNA_SMARTer_exprs.npz']
            self.rna_labels = ['data_MOp/YaoEtAl_RNA_snRNA_10X_v3_B_cellTypes.txt',\
                                'data_MOp/YaoEtAl_RNA_snRNA_10X_v3_A_cellTypes.txt',\
                                'data_MOp/YaoEtAl_RNA_snRNA_10X_v2_cellTypes.txt',\
                                'data_MOp/YaoEtAl_RNA_snRNA_SMARTer_cellTypes.txt',\
                                'data_MOp/YaoEtAl_RNA_scRNA_10X_v3_cellTypes.txt',\
                                'data_MOp/YaoEtAl_RNA_scRNA_10X_v2_cellTypes.txt',\
                                'data_MOp/YaoEtAl_RNA_scRNA_SMARTer_cellTypes.txt']
            self.atac_paths = ['data_MOp/YaoEtAl_ATAC_exprs.npz',\
                                'data_MOp/YaoEtAl_snmC_exprs.npz']
            self.atac_labels = ['data_MOp/YaoEtAl_ATAC_cellTypes.txt',\
                                'data_MOp/YaoEtAl_snmC_cellTypes.txt']
            self.rna_protein_paths = []
            self.atac_protein_paths = []
            
            # Training config            
            self.batch_size = 256
            self.lr_stage1 = 0.01
            self.lr_stage3 = 0.001
            self.lr_decay_epoch = 20
            self.epochs_stage1 = 10
            self.epochs_stage3 = 10
            self.p = 0.8
            self.embedding_size = 64
            self.momentum = 0.9
            self.center_weight = 20
            self.with_crossentorpy = True
            self.seed = 1
            self.checkpoint = '' 
            
        elif DB == "db4_control":
            self.number_of_class = 7 # Number of cell types in CITE-seq data
            self.input_size = 17668 # Number of common genes and proteins between CITE-seq data and ASAP-seq
            self.rna_paths = ['data/citeseq_control_rna.npz'] # RNA gene expression from CITE-seq data
            self.rna_labels = ['data/citeseq_control_cellTypes.txt'] # CITE-seq data cell type labels (coverted to numeric) 
            self.atac_paths = ['data/asapseq_control_atac.npz'] # ATAC gene activity matrix from ASAP-seq data
            self.atac_labels = ['data/asapseq_control_cellTypes.txt'] # ASAP-seq data cell type labels (coverted to numeric) 
            self.rna_protein_paths = ['data/citeseq_control_adt.npz'] # Protein expression from CITE-seq data
            self.atac_protein_paths = ['data/asapseq_control_adt.npz'] # Protein expression from ASAP-seq data
            
            # Training config            
            self.batch_size = 256
            self.lr_stage1 = 0.01
            self.lr_stage3 = 0.01
            self.lr_decay_epoch = 20
            self.epochs_stage1 = 40
            self.epochs_stage3 = 20
            self.p = 0.8
            self.embedding_size = 64
            self.momentum = 0.9
            self.center_weight = 1
            self.with_crossentorpy = True
            self.seed = 1
            self.checkpoint = '' 

        elif DB == "CoupleCocData":
            self.number_of_class = 3 # Number of cell types in CITE-seq data
            self.input_size = 1000 # Number of common genes and proteins between CITE-seq data and ASAP-seq
            self.rna_paths = ['coupleCocData/example1A_S_label.npz'] # RNA gene expression from CITE-seq data
            self.rna_labels = ['coupleCocData/example1A_S.txt'] # CITE-seq data cell type labels (coverted to numeric) 
            self.atac_paths = ['coupleCocData/example1A_T_label.npz'] # ATAC gene activity matrix from ASAP-seq data
            self.atac_labels = ['coupleCocData/example1A_T.txt'] # ASAP-seq data cell type labels (coverted to numeric) 
            self.rna_protein_paths = [] # Protein expression from CITE-seq data
            self.atac_protein_paths = [] # Protein expression from ASAP-seq data
            
            # Training config            
            self.batch_size = 256
            self.lr_stage1 = 0.01
            self.lr_stage3 = 0.01
            self.lr_decay_epoch = 20
            self.epochs_stage1 = 40
            self.epochs_stage3 = 20
            self.p = 0.8
            self.embedding_size = 64
            self.momentum = 0.9
            self.center_weight = 1
            self.with_crossentorpy = True
            self.seed = 1
            self.checkpoint = '' 



        elif DB == "HematopoiesisData":
            self.number_of_class = 26 # Number of cell types in CITE-seq data
            self.input_size = 15714 # Number of common genes and proteins between CITE-seq data and ASAP-seq
            self.rna_paths = ['HematopoiesisData/rna_filtered_t.npz'] # RNA gene expression from CITE-seq data
            self.rna_labels = ['HematopoiesisData/rna_labels.txt'] # CITE-seq data cell type labels (coverted to numeric) 
            self.atac_paths = ['HematopoiesisData/atac_filtered_t.npz'] # ATAC gene activity matrix from ASAP-seq data
            self.atac_labels = ['HematopoiesisData/atac_labels.txt'] # ASAP-seq data cell type labels (coverted to numeric) 
            self.rna_protein_paths = [] # Protein expression from CITE-seq data
            self.atac_protein_paths = [] # Protein expression from ASAP-seq data
            
            # Training config            
            self.batch_size = 256
            self.lr_stage1 = 0.01
            self.lr_stage3 = 0.01
            self.lr_decay_epoch = 20
            self.epochs_stage1 = 40
            self.epochs_stage3 = 20
            self.p = 0.8
            self.embedding_size = 64
            self.momentum = 0.9
            self.center_weight = 1
            self.with_crossentorpy = True
            self.seed = 1
            self.checkpoint = '' 


            

        



