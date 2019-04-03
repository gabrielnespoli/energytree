# Create the dataframe used in the regression/classification using graphs as features
# The graphs (x-data) are contained in [graph_folder] and the Y-data is present in [y.file]'

generate_graph_data <- function(graph_folder='./data/NKI_Rockland/',
                                y.filename='./data/clinical_information.txt',
                                y.column = 'WASI_FULL_4',
                                output_filename = 'data.RData',
                                output_folder = "./data/"
                                )
  {
  
  require(igraph)
  y.file = read.csv(y.filename)[c('Subject',y.column)]
  all.und.graphs <- list()
  all.dir.graphs <- list()
  all.ys <- list()
  i = 1
  graph.files = list.files(graph_folder, pattern='_DTI_connectivity_matrix_file.txt', full.names = FALSE) # retrieve all the filenames that keep the graph data
  for(filename in graph.files) {
    id = strsplit(filename, split = '_', fixed = TRUE)[[1]][1] # get the patient id from the filename
    
    und.filename = paste(graph_folder,filename, sep='')
    undirected = as.matrix(read.table(und.filename)) # 'paste' concatenates folder and filename. Read the undirected graph as a matrix
    undirected = graph_from_adjacency_matrix(adjmatrix = undirected,mode = 'undirected', diag=FALSE, weighted=TRUE)
    
    result = tryCatch({
      dir.filename = paste(graph_folder,id, '_fcMRI_noGSR_connectivity_matrix_file.txt', sep='')
      directed = as.matrix(read.table(dir.filename)) # Read the directed graph as a matrix
      directed = graph_from_adjacency_matrix(adjmatrix = directed,mode = 'directed', diag=FALSE, weighted=TRUE)
    }, error = function(e) {
      dir.filename = paste(graph_folder,id, '_fcMRI_GSR_connectivity_matrix_file.txt', sep='')
      directed = as.matrix(read.table(dir.filename)) # Read the directed graph as a matrix
      directed = graph_from_adjacency_matrix(adjmatrix = directed,mode = 'directed', diag=FALSE, weighted=TRUE)
    }
    )
    
    y = y.file[y.file$Subject == as.numeric(id), y.column]
    
    if(!is.na(y)) {
      all.und.graphs[[i]] = undirected
      all.dir.graphs[[i]] = directed
      all.ys[[i]] = y
      i = i + 1
    }
    
  }
  
  # creates an empty list that will store the data
  # x will store the graph data (the covariate) and y the response variable 
  data <- list()
  data[['Y']] = as.numeric(all.ys)
  data[['V1']] = all.und.graphs
  data[['V2']] = all.dir.graphs
  
  save(data, file = paste( output_folder, output_filename , sep = ''))
  
}

# generate_graph_data()

generate_graph_data(graph_folder='../data/NKI_Rockland/',
                    y.filename='../data/clinical_information.txt',
                    y.column = 'TFEQ_SUSC_HUNGER',
                    output_filename = 'data2.RData',
                    output_folder = "../data/"
)

load('../data/data2.RData')
