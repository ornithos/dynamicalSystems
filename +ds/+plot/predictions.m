function predictions(obj, pwindows, horizOverlay, varargin)
% plot.predictions(obj, pwindows, horizOverlay, titlePrefix)
% plot prediction of dynamicalSystem object
%
%   INPUTS:
%      obj           - dynamicalSystems object
%      pwindows      - (matrix) prediction windows (list of, ie predict k-ahead
%                       for k in pwindows).
%      horizOverlay  - (matrix) background colour of plot: same length as T
%      titlePrefix   - (str) prefix of title of each plot (eg. 'system 3:')

      % Setup
      assert(isa(obj, 'ds.dynamicalSystem'), 'input object not a dynamicalSystems object');
      isBatch       = isa(obj, 'ds.dynamicalSystemBatch');
      
      if nargin < 2 || isempty(pwindows)
          pwindows = 1;
      end
      assert(all(arrayfun(@(x) utils.is.scalarint(x,1), pwindows)), 'pwindows must be positive scalar int (matrix)');
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

      optsDefault.titlePrefix = '';
      if isBatch
        optsDefault.nn        = 1:obj.d.n;
      else
        optsDefault.nn        = 1;
      end
      optsDefault.verbose     = true;
      optsDefault.figure      = [];
      opts           = utils.base.processVarargin(varargin, optsDefault);
          
      
          
      cols          = zeros(obj.d.y,3);
      cnums         = [2 4 5 3 6 7 1];
      for kk = 1:obj.d.y; cols(kk,:) = utils.plot.varyColor2(cnums(kk)); end
      litecols      = utils.plot.colortint(cols, 0.5);

      %% Plot
      if isempty(opts.figure)
          figure;
      else
          figure(opts.figure);
      end
      
      nPreds        = numel(pwindows);
      spdims        = utils.plot.subplotdims(nPreds);

      % get data
      if isBatch && opts.verbose; fprintf('Calculating: fitted ys...'); end
      [impy, impPy] = obj.impute_y('variance', true, 'smooth', true, 'nRng', opts.nn);
      predvals      = cell(nPreds,1);
      
      if isBatch && opts.verbose; fprintf('\b\b\b, predicted values, ...'); end
      for jj = 1:nPreds
          predvals{jj} = obj.getPredictedValues(pwindows(jj), [], opts.nn);
          if ~isBatch; predvals{jj} = {predvals{jj}}; end
      end
      if isBatch && opts.verbose; fprintf('\b\b\b\b\b. Done.\n'); end
      
      % data may be cell if dsBatch
      if isBatch
          y       = obj.y;
          d       = obj.ambientDimension;
      else
          opts.nn = 1;
          impy = {impy}; impPy = {impPy}; y = {obj.y};
          horizOverlay = {horizOverlay};
          d       = obj.d.y;
      end
      
      for nn = opts.nn
          for jj = 1:nPreds
              subplot(spdims(1), spdims(2), jj);

              for kk = 1:d(nn)
    %                       plot(obj.y(jj,:)', ':', 'Color', litecols(jj,:)); hold on;
                  plot(impy{nn}(kk,:)', '-', 'Color', litecols(kk,:)); hold on;
              end
              stdy          = sqrt(cell2mat(cellfun(@(x) diag(x)', impPy{nn}, 'Un', 0)))';

              for kk = 1:d(nn)
                  nanidx          = find(isnan(y{nn}(kk,:)));
                  if isempty(nanidx); continue; end
                  nanminmax       = nanidx(1):nanidx(end);
                  utils.plot.confidenceInterval(nanminmax, impy{nn}(kk,nanminmax), stdy(kk,nanminmax), [], 'facecolor', cols(kk,:));
              end
              
              ix       = (obj.d.T(nn) - size(predvals{jj}{nn},2)+1):obj.d.T(nn);
              smse     = nanmean((y{nn}(:,ix) - predvals{jj}{nn}).^2,2)./nanvar(y{nn}, [], 2);   % smse by channel (STNDD BY VAR(*all* Y))
              for kk = 1:d(nn)
                  plot(ix, predvals{jj}{nn}(kk,:), '-', 'Color', cols(kk,:)); 
              end
              hold off;

              if ~isempty(horizOverlay)
                  cmap           = flipud([1 0.86 0.91; 0.86 0.91 1]);
                  if size( horizOverlay{nn}, 2) == 1;  horizOverlay{nn} =  horizOverlay{nn}'; end
                  utils.plot.dataShadeVertical(1:(obj.d.T(nn)), horizOverlay{nn}, cmap, 'edgedetect', 'manual', 'edgemanual', 3);
              end
              
              if numel(opts.nn) > 1
                  idtxt = sprintf('Series %d', nn);
              else
                  idtxt = '';
              end
              title(sprintf('%s%s (%d-step). Average SMSE %.2f%%', opts.titlePrefix, idtxt, pwindows(jj), nanmean(smse)*100));
          end
          
          if nn < opts.nn(end)
              pause;
          end
      end

end
