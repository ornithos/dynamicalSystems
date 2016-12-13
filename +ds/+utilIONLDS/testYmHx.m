pb = utils.base.objProgressBar;
pb.newProgressBar('progress: ', 30, 'showElapsed', true);

ymHxMC  = zeros(dsINL.d.T, dsINL.d.y);
outerMC = zeros(dsINL.d.y, dsINL.d.y, dsINL.d.T);
for tt = 1:dsINL.d.T
    ns = 20000;
    emiParams = dsINL.par.emiNLParams;
    
    cholS   = chol(dsINL.infer.smooth.sigma{tt});
    ymHx    = zeros(dsINL.d.y, 1);
    outer   = zeros(dsINL.d.y,dsINL.d.y);
    for ii = 1:ns
        x      = dsINL.infer.smooth.mu(:,tt) + cholS * randn(dsINL.d.x,1);
        hx     = ds.utilIONLDS.gen_sigmoid(emiParams.C * x, emiParams.eta);
        y      = dsINL.y(:,tt);
        df     = y - hx;
        ymHx   = ((ii-1)*ymHx + df)./ii;
        outer  = ((ii-1)*outer + df*df')./ii;
    end
    
    pb.print(tt/dsINL.d.T);
    ymHxMC(tt,:) = ymHx;
    outerMC(:,:,tt) = outer;
end

pb.finish;