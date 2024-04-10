function pred = classf(X1train,ytrain,X1test)
Xtrain = table(X1train,ytrain,'VariableNames',{'Speed','Drag'});
Xtest = table(X1test,'VariableNames',{'Speed'});
modelspec = 'Drag ~ Speed';
mdl = fitglm(Xtrain,modelspec,'Distribution','normal');
yfit = predict(mdl,Xtest);
pred = (yfit > 0.5);
end