function [nx, ny, nz, matrix3D] = loadmat3d(filename)
  % Open the file
  fid = fopen(filename, 'r');
  
  if fid == -1
      error('Failed to open file: %s', filename);
  end
  
  % Read the dimensions of the matrix
  dimensions = fscanf(fid, '%d %d %d', [1 3]);
  nx = dimensions(1);
  ny = dimensions(2);
  nz = dimensions(3);
  
  % Initialize the 3D matrix
  matrix3D = zeros(nz, ny, nx);
  
  % Read the matrix data
  for i = 1:nx
      fgetl(fid);
      for j = 1:ny
          line = fgetl(fid);
          data = sscanf(line, '%f');
          matrix3D(:,j,i) = data;
      end
  end

  matrix3D = permute(matrix3D, [3 2 1]);
  % for j = 1:ny
  %     % fgetl(fid);
  %     for k = 1:nz
  %         % Read one line of data
  %         line = fgetl(fid);
  %         % fprintf(line);
  %         line
  %         % Convert the line to numerical data
  %         % data = sscanf(line, '%f');
  %         % Store the data in the 3D matrix
  %         % matrix3D(:, j, k) = data;
  %     end
  %     % Skip the empty line between layers
  %     % fgetl(fid);
  % end
  
  % Close the file
  fclose(fid);
end
