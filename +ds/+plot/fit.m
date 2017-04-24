function fit(obj, pforward, horizOverlay)
% plot.fit(obj, pforward, horizOverlay)
% plot fit of dynamicalSystem object
%
%   * Fitted values
%   * Includes credible intervals of inferred outputs (where NaN)
%   * shows prediction of time series after end (# = pforward)
%   * Can overlay some 1d matrix of values to colour background
%   (horizOverlay).

      % Setup
      assert(isa(obj, 'ds.dynamicalSystem'), 'input object not a dynamicalSystems object');
      if nargin < 2 || isempty(pforward)
          pforward = 1;
      end
      assert(utils.is.scalarint(pforward, 1), 'pforward must be a positive scalar int');
      if nargin < 3
          horizOverlay = [];
      else
          if iscell(horizOverlay)
              assert(all(cellfun(@numel, horizOverlay) == obj.d.T), 'horizOverlay wrong size for (some) time series');
          else
              assert(~isa(obj, 'ds.dynamicalSystemBatch'), 'dynamicalSystemBatch requires a cell input of horizOverlay');
              assert(numel(horizOverlay) == obj.d.T, 'horizOverlay wrong size for tinme series');
          end
      end

      cols          = zeros(obj.d.y,3);
      cnums         = [2 4 5 3 6 7 1];
      for kk = 1:obj.d.y; cols(kk,:) = utils.plot.varyColor2(cnums(kk)); end
      litecols      = utils.plot.colortint(cols, 0.8);

      %% Plot
      figure
      
      % get data
      [impy, impPy] = obj.impute_y('variance', true, 'smooth', true);
      futurY        = obj.getPredictFreeRun(obj.d.T, pforward);
      
      % data may be cell if dsBatch
      if isa(obj, 'ds.dynamicalSystemBatch')
          Nmodels = obj.d.n;
          y       = obj.y;
      else
          Nmodels = 1;
          impy = {impy}; impPy = {impPy}; futurY = {futurY}; y = {obj.y};
          horizOverlay = {horizOverlay};
      end
      
      for nn = 1:Nmodels
          for jj = 1:obj.d.y
    %               plot(obj.y(jj,:)', ':', 'Color', litecols(jj,:)); hold on;
              plot(impy{nn}(jj,:)', '-', 'Color', cols(jj,:)); hold on;
              plot(1:obj.d.T(nn)+pforward, futurY{nn}(jj,:), 'Color', litecols(jj,:));
          end
          hold off;
          stdy          = sqrt(cell2mat(cellfun(@(x) diag(x)', impPy{nn}, 'Un', 0)))';
          for jj = 1:obj.d.y
              nanidx          = find(isnan(y{nn}(jj,:)));
              if isempty(nanidx); continue; end
              nanminmax       = nanidx(1):nanidx(end);
              utils.plot.confidenceInterval(nanminmax, impy{nn}(jj,nanminmax), stdy(jj,nanminmax), [], 'facecolor', cols(jj,:));
          end
          if ~isempty(horizOverlay)
              cmap           = flipud([1 0.86 0.91; 0.86 0.91 1]);
              if size( horizOverlay{nn}, 2) == 1;  horizOverlay{nn} =  horizOverlay{nn}'; end
              utils.plot.dataShadeVertical(1:(obj.d.T(nn)), horizOverlay{nn}, cmap, 'edgedetect', 'manual', 'edgemanual', 3);
          end
          
          if nn < Nmodels
              pause;
          end
          
          if nModels > 1
              title('Series %d', nn);
          end
      end
end
