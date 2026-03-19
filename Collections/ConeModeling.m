function model=ConeModeling(model,alphaPos,LinConst,ConeIndices)

% model.modelsense='min';
% model.obj=full(sparse(1,alphaPos,-1,1,size(LinConst,2))); % Objective

model.modelsense='max';
model.obj=full(sparse(1,alphaPos,1,1,size(LinConst,2))); % Objective

model.A=LinConst; % Equality
model.rhs=zeros(size(LinConst,1),1);
model.sense='=';

% Inequality constraints ()
for x=1:numel(ConeIndices)
    model.cones(x).index=ConeIndices{x};
end

end

