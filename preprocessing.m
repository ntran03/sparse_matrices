matrixName = "";
%allMatrix(matrixName)
for n = 0:168
    strcat("epb3_", num2str(n))
    csr_matrix(strcat("epb3_", num2str(n)))
end
%% split one large matrix and rank each chunk
%requires individual determination of how many to split the original matrix
%into based on matrix size relative to 8gb
function all = allMatrix(matrixName)
    %split matrix into 16ths
    matrix = openMatrixFile(matrixName);
    mSize = size(matrix);
    x = mSize(1);
    %puts the 16ths into a folder under the name of the matrix
    splitMegaMatrix(openMatrixFile(matrixName),x/16, matrixName(1:end-4))
    %run main on each of the matrices inside the folder
    numMatrices = 256
    for k = 0:numMatrices-1
        tempName = matrixName(1:end-4)+"_" + num2str(k) + ".mat";
        
        main(openMatrixSubFile(tempName,matrixName(1:end-4)), erase(tempName, ".mat"))
    end
end
%% split megamatrix
function numMatrices = splitMegaMatrix(matrix, n, folder)
    mkdir(folder)
    matrixID = 0
    mSize = size(matrix)
    x = mSize(1)
    y = mSize(2)
    C = mat2cell(matrix, ones(1,x/n)*n, ones(1,y/n)*n);
    %celldisp(C)
    %for each cell in C, save it as a .mat file with name
    %matrix(matrixID).mat
    for i = 1:(x/n)
        for j = 1:(y/n)
            tempMatrix = C(i,j);
            save([folder+"/"+folder+"_"+num2str(matrixID) + ".mat"], 'tempMatrix')
            matrixID = matrixID + 1
        end
    end
    numMatrices = x*y/(n * n);
end
%% splitMatrix (done)
%splits matrix into nxn matrices, reading like a books
%saves all matrices as .mat files to folder
%returns numMatrices, the number of matrices made
function numMatrices = splitMatrix(matrix, n, folder)
    folder
    mkdir(folder)
    matrixID = 0
    mSize = size(matrix)
    x = mSize(1)
    y = mSize(2)
    C = mat2cell(matrix, ones(1,x/n)*n, ones(1,y/n)*n);
    %celldisp(C)
    %for each cell in C, save it as a .mat file with name
    %matrix(matrixID).mat
    for i = 1:(x/n)
        for j = 1:(y/n)
            tempMatrix = C(i,j);
            save([folder+"/"+"matrix"+"_"+num2str(matrixID) + ".mat"], 'tempMatrix')
            matrixID = matrixID + 1
        end
    end
    numMatrices = x*y/(n * n);
end
%% padding (done)
%add rows and columns of 0s as necessary to get to row/column being
%multiples of n, the goal matrix dimension
function X = padding(matrix, n)
    %pad the array
    mSize = size(matrix);
    x = mod(mSize(1), n);
    if  x == 0
        
    else
        x = n-mod(mSize(1),n);
    end
    y = mod(mSize(2), n)
    if  y == 0
        
    else
        y = n-mod(mSize(2),n);
    end
    X = padarray(matrix, [x y], 0, "post");
end

%% diagonal/ranking (done)
%define the band of diagonalness
%find the number of nonzeroes in the band
%divide by total number of nonzeroes
%output diagonalness
%rank 1, rank 4, rank 16: output 3 values for each matrix unit
%one for each
%BIG EDIT: multiply each ranking by the # of nonzero values in the block
function rankings = rankMatrix(matrix, band, folder) 
    %load matrix
    currentMatrix = load(fullfile(folder,matrix), 'tempMatrix');
    current = cell2mat(getfield(currentMatrix, 'tempMatrix'));
    %retrieve matrix ID
    pat = digitsPattern;
    matrixID = extract(matrix,pat);
    %get rank_1
    matrix1 = current;
    %rank the chunk
    temp = rankChunk(matrix1, 256, band);
    rank_1 = temp*100;
    %get rank_4
    ranklist_4 = [];
    %split matrix into 4
    matrix4 = mat2cell(current, [128 128], [128 128]);
    for i = 1:2
        for j = 1:2
            %rank the chunk and add to the rank_4 array
            temp = rankChunk(cell2mat(matrix4(i,j)), 128, band);
            ranklist_4 = [ranklist_4, temp];
        end
    end
    rank_4 = mean(ranklist_4)*100;
    %get rank_16
    ranklist_16 = [];
    matrix16 = mat2cell(current, [64 64 64 64], [64 64 64 64]);
    for i = 1:4
        for j = 1:4
            %rank the chunk and add to the rank_16 array
            temp = rankChunk(cell2mat(matrix16(i,j)), 64, band);
            ranklist_16 = [ranklist_16, temp];
        end
    end
    rank_16 = mean(ranklist_16)*100;
    %final array to return
    rankings = [matrixID rank_1*nnz(matrix1) rank_4*nnz(matrix1) rank_16*nnz(matrix1) nnz(matrix1)];
end

%returns the ranking for a given piece of the matrix
%band width = w, dimension of matrix = n
function temprank = rankChunk(matrix, n, w)
    nonzero = nnz(matrix);
    diagNonzero = 0;
    %for a standard width
    for i = -(w/2):(w/2)
        diagNonzero = diagNonzero + nnz(diag(matrix, i));
    end
    if nonzero ~= 0
        temprank = diagNonzero/nonzero;
    else
        temprank = 0;
    end
end

%% open files (done)
%read the downloaded matrix
%takes in string that is the matrix name
function A = openMatrixFile(name)
    %in theory, opens and reads the matrix
    fullStructure = load(name, 'Problem');
    T = getfield(fullStructure, 'Problem')
    A = getfield(T,'A')
    
end
%% open chunked matrix files
%when I had to break down the matrices into more managable chunks
function T = openMatrixSubFile(name, folder)
    fileName = fullfile(folder, name);
    matrix = load(fileName, 'tempMatrix');
    T = getfield(matrix, 'tempMatrix');
    T = cell2mat(T);
end
%% csr format
function csr_matrix(filen)
    dirr = strcat("/Users/nicol/Documents/sparsity/epb3/",filen)
    %file_dir = dir('/Users/nicol/Documents/sparsity/2cubes_sphere/2cubes_sphere_2/*.mat');
    file_dir = dir(strcat('/Users/nicol/Documents/sparsity/epb3/',filen,'/*.mat'));
    read_dir = dirr;
    dir_name = dirr;
    foldername = dirr;
%the way this file works is it breaks down a whole matrix, but what I need
%to get out of it is just the csr portion- look for where you get offsets,
%indices, etc, and then use my own code to iterate through the matrices I
%have and save them in a separate folder.
    
        %okay so what happens is each big matrix is broken into blocks, and
        %then each block gets a line in each file
        %however, I kind of want one csr thing per each matrix since I
        %broke mine down already.
        %how I can do this is first open one file and use it to set the
        %block width, and iterate through each of the things in the folder
    CSR_offsets = fopen(strcat(dir_name,'/',filen,'_CSR_OFFSETS'),'w');
    CSR_col_indices = fopen(strcat(dir_name,'/',filen,'_CSR_COL_INDICES'),'w');
    CSR_values = fopen(strcat(dir_name,'/',filen,'_CSR_VALUES'),'w');
    
    for k=1:length(file_dir)
        file_name = file_dir(k).name;
        file_name

        % ------ output (I want them to go in three files that will sit in
        % the original folder with the uncompressed matricessss

        % 2. CSR

       
        %make a separate thing to load my files!
        %sike ust use the openmatrixsubfile thing
        
        %Data = load(strcat(read_dir,file_name));
        %A = Data.Problem.A;
        A = openMatrixSubFile(file_name, foldername);
        A_size = size(A); % a 2x2 matrix including the dimension of A
        % create CoperSim input files by compressing a
        % block (C{i,j}) into different formats. one block
        % will be writen as a single line in the correspoinding file(s)


        %THIS IS THE ACTUAL CSR STUFF- BY BLOCK THEY MEAN
        %WHAT THEY CREATED, BUT I WILL DIRECTLY USE THE MATRIX I HAVE
        %I wonder if it will work without perfect sized things
        % 2. CSR
      
        offset = 0;
        block_width = A_size(1);
        for ii = 1:block_width
            offset = offset + nnz(A(ii,:));
            %why are there two things in nnz
            %oh it's because actually I do not know haha
            fprintf(CSR_offsets,'%d ', offset); % number of nz values in a row
            for jj = 1:block_width
                if (A(ii,jj) ~= 0) % write all the nz values and corresponding col indices per each row
                    fprintf(CSR_col_indices,'%d ', jj); % jj is the col index of a value
                    fprintf(CSR_values,'%d ', 1); % C{i,j}(ii,jj)); % writie the value itself
                end
            end
         end
         fprintf(CSR_offsets,'\n');
         fprintf(CSR_col_indices,'\n');
         fprintf(CSR_values,'\n');
    end
         fclose(CSR_offsets);
         fclose(CSR_col_indices);
         fclose(CSR_values);
 disp('done')
 end
                       
        
        


       
%% main
%opens and ranks all matrices
%end goal is to have something where the matrices and their corresponding
%rankings are stored? ask in meeting
function allRankings = main(matrix, name)
    %split matrix and get # of matrices
    matrix = padding(matrix, 256);
    numMatrices = splitMatrix(matrix, 256, name);
    rankName = name+'/rankings_'+name+'.txt'
    fid = fopen(rankName,'wt');
    %initialize an array for the rankings (by 5 because 5 things in output
    %vector
    allRankings = zeros(numMatrices,5);
    for k = 0:numMatrices-1
        matrixName = "matrix_" + num2str(k) + ".mat";
        %rank the split matrices
        ranking = rankMatrix(matrixName, 8, name);
        %add to the rankings array
        allRankings(k+1,:) = ranking;
        for i = 1:4
            fprintf(fid, '%s ',ranking(i));
        end
        fprintf(fid, '%s\n', ranking(5));
    end    
    fclose(fid);
end

