
% initial timeseries volume info

function fsol = InitFlowHDF5(dim)
  fsol.gname = '/';
  fsol.sname = 'time';
  fsol.rank = dim;
  fsol.count = ones(1,dim);
  fsol.stride = ones(1,dim);
  fsol.offset = zeros(1,dim);
end
