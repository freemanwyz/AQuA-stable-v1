function res = evt2lmkProp1(datS,lmkMsk)
% distances and directions between events and landmarks
% Multiple threshold frontier based

[H,W,T] = size(datS);
nLmk = numel(lmkMsk);

thrRg = 0.2:0.1:0.8;
% thrRg = 0.5;
chgToward = zeros(1,nLmk);
chgAway = zeros(1,nLmk);
chgTowardBefReach = zeros(1,nLmk);
chgAwayAftReach = zeros(1,nLmk);
for kk=1:numel(thrRg)
%     fprintf('%d\n',kk)
    evt0 = datS>thrRg(kk);
    loc0 = find(evt0>0);
    
    % should not happen
    if isempty(loc0)
        continue
    end
    
    % impute missed frames
    % some times some frame could be missed due to some post processing
    [~,~,it] = ind2sub([H,W,T],loc0);
    tRg = min(it):max(it);
    if numel(tRg)==1
        continue
    end    
    for tt=1:numel(tRg)
        ixSel = it==tRg(tt);
        if isempty(ixSel)
            evt0(:,:,tRg(tt)) = evt0(:,:,tRg(tt)-1);
        end
    end
    
    % distance to landmarks
    % use the center of the landmark
    % use geodesic distance
    % % if landmark outside event, use a pixel in event closest to landmark as landmark
    % if landmark outside event, use Euclidean distances
    D = cell(nLmk,1);
    evt0s = squeeze(sum(evt0,3)>0);
    for ii=1:nLmk
        msk00 = lmkMsk{ii};
        [h0,w0] = find(msk00>0);
        h00 = mean(h0); w00 = mean(w0);
        msk00 = zeros(H,W); msk00(max(round(h00),1),max(round(w00),1)) = 1;
        if sum(evt0s(msk00>0))==0
            [h1,w1] = find(evt0s>0);
            dist00 = sqrt((h1-h00).^2+(w1-w00).^2);
            %[~,ix] = min(dist00);
            %msk00 = zeros(H,W); msk00(h1(ix),w1(ix)) = 1;
            tmp = zeros(H,W);
            tmp(evt0s>0) = dist00;
        else
            tmp = bwdistgeodesic(evt0s,msk00>0);
        end        
        D{ii} = tmp;
    end
    
    % frontier tracking
    bdLst = cell(numel(tRg),1);
    lblMap = zeros(H,W,numel(tRg));    
    ccLst = cell(numel(tRg),1);
    for ii=1:numel(tRg)
        xCur = evt0(:,:,tRg(ii));
        [B,L] = bwboundaries(xCur,8,'noholes');
        bdLst{ii} = B;
        lblMap(:,:,ii) = L;
        ccLst{ii} = label2idx(L);
    end

    dxAllPos = zeros(numel(tRg),nLmk);
    dxAllNeg = zeros(numel(tRg),nLmk);
    tReach = nan(1,nLmk);
    for ii=2:numel(tRg)
        lblCur = lblMap(:,:,ii);
        for jj=1:nLmk
            % !! suffient coverage of landmark
            %tmp = lmkMsk{jj}; 
            %n00 = sum(tmp(:)>0);
            n11 = sum(lblCur(lmkMsk{jj}>0));
            insideLmk = n11>0;
            %insideLmk = n11>n00*0.5;
            %insideLmk = n11>10 || n11>n00*0.3;
            if insideLmk && isnan(tReach(jj))
                tReach(jj) = ii;
            end
        end
        
        ccCur = ccLst{ii};
        lblPre = lblMap(:,:,ii-1);
        for jj=1:numel(ccCur)
            % cc in previous frame that connect to this cc
            cc0 = ccCur{jj};
            lblSel = unique(lblPre(cc0));
            lblSel = lblSel(lblSel>0);
            if isempty(lblSel)
                continue
            end
            
            % previous boundary, the starting point of propagation
            % we use some ad hocs:
            % smaller pre area has higher distance penalty
            % multiple pre area could compete
            % if a previous cc is too small itself and/or relative to current cc, ignore it
            % may not be the correct one, but looks more comfortable
            bdPre = [];
            bdPreWt = [];
            n0c = numel(cc0);
            for uu=1:numel(lblSel)
                n0 = numel(ccLst{ii-1}{lblSel(uu)});
                if n0>n0c/5
                    tmp = bdLst{ii-1}{lblSel(uu)};
                    tmp = sub2ind([H,W],tmp(:,1),tmp(:,2));
                    bdPre = [bdPre;tmp]; %#ok<AGROW>
                    bdPreWt = [bdPreWt;ones(numel(tmp),1)/n0]; %#ok<AGROW>
                end
            end
            if isempty(bdPre)
                continue
            end
            
            % current boundary, the ending point of propagation
            % or use all increased pixels?
            % do not include boundary that is active in previous frame
            % we only use the incresing signals

            bdCur = ccLst{ii}{jj};            
            bdCur = bdCur(lblPre(bdCur)==0);            
            if isempty(bdCur)
                continue
            end
            
            % link each pixel in bdCur to a pixel in bdPre
            % for each landmark, find the distance change for each pair
            % positive change is away from landmark
            % if pixCur contains landmark, it is treated as two parts
            dxPos = zeros(numel(bdCur),nLmk);
            dxNeg = zeros(numel(bdCur),nLmk);
            [h0,w0] = ind2sub([H,W],bdCur);
            [h1,w1] = ind2sub([H,W],bdPre);
            for uu=1:numel(bdCur)
                % closest starting frontier point                                
                d00 = sqrt((h0(uu)-h1).^2+(w0(uu)-w1).^2);
                d01 = d00.*bdPreWt;
                [~,ix] = min(d01);
                d00min = d00(ix);                
                
                % find path between points
                h0a = h0(uu); w0a = w0(uu); h1a = h1(ix); w1a = w1(ix);                
                wGap = (w1a-w0a)/max(round(d00min),1);
                hGap = (h1a-h0a)/max(round(d00min),1);
                hx = round(h0a:hGap:h1a);
                wx = round(w0a:wGap:w1a);
                if h0a==h1a && w0a==w1a
                    hx = h0a; wx = w0a;
                elseif h0a==h1a
                    hx = ones(1,numel(wx))*h0a;
                elseif w0a==w1a
                    wx = ones(1,numel(hx))*w0a;
                end                
                hwx = sub2ind([H,W],hx,wx);

                % propagation distance w.r.t landmarks
                for vv=1:nLmk
                    D0 = D{vv};
                    dp0 = D0(hwx);
                    dp0Min = min(dp0);
                    dxPos(uu,vv) = max(D0(h1a,w1a)-dp0Min,0);  % toward
                    dxNeg(uu,vv) = max(D0(h0a,w0a)-dp0Min,0);  % away               
                end
            end
            
            if 0
                lmkSel = 2;
                tmp1 = zeros(H,W); bd1 = bdCur(dxPos(:,lmkSel)>0); tmp1(bd1) = 1;
                tmp2 = zeros(H,W); bd2 = bdCur(dxNeg(:,lmkSel)>0); tmp2(bd2) = 1;
                tmp3 = lmkMsk{vv}*0.3; tmp3(bdPre) = 1;
                tmp = cat(3,tmp1,tmp2,tmp3); figure;imshow(tmp)                
                text(20,20,sprintf('Toward %f - Away %f',sum(dxPos(:,lmkSel)),sum(dxNeg(:,lmkSel))),'Color','y');
                %pause(2); 
                keyboard
                close
            end
            
            dxAllPos(ii,:) = dxAllPos(ii,:) + sum(dxPos,1);
            dxAllNeg(ii,:) = dxAllNeg(ii,:) + sum(dxNeg,1);
        end
    end
    chgToward = chgToward + sum(dxAllPos,1);
    chgAway = chgAway + sum(dxAllNeg,1);
    for ii=1:nLmk
        if ~isnan(tReach(ii))
            t1 = min(tReach(ii)+1,numel(tRg));
            chgTowardBefReach(ii) = chgTowardBefReach(ii) + sum(dxAllPos(1:tReach(ii),ii));
            chgAwayAftReach(ii) = chgAwayAftReach(ii) + sum(dxAllNeg(t1:end,ii));
        end
    end
end

res.chgToward = chgToward;
res.chgAway = chgAway;
res.chgTowardBefReach = chgTowardBefReach;
res.chgAwayAftReach = chgAwayAftReach;

end






