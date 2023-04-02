function [Uhis,Fhis] = PathAnalysis(truss,angles,AnalyInputOpt)
    tol = 1e-6; MaxIter = 50; 
    Node = truss.Node;
    if ~isfield(truss,'U0'), truss.U0 = zeros(3*size(truss.Node,1),1); end  
    U = truss.U0;
    
    if strcmpi(AnalyInputOpt.LoadType, 'Force')
        MaxIcr = AnalyInputOpt.MaxIcr;
        b_lambda = AnalyInputOpt.InitialLoadFactor;
        Uhis = zeros(3*size(Node,1),MaxIcr);
        FreeDofs = setdiff(1:3*size(Node,1),truss.FixedDofs);
        lmd = 0; icrm = 0; MUL = [U,U];
        Fhis = zeros(MaxIcr,1);
        if isfield(AnalyInputOpt,'Load')
            F = AnalyInputOpt.Load;
        end
        while icrm<MaxIcr && ~AnalyInputOpt.StopCriterion(Node,U,icrm)
            icrm = icrm+1;
            iter = 0; err = 1;
            fprintf('icrm = %d, lambda = %6.4f\n',icrm,lmd);
            if isfield(AnalyInputOpt,'AdaptiveLoad')
                F = AnalyInputOpt.AdaptiveLoad(Node,U,icrm);
            end
            while err>tol && iter<MaxIter
                iter = iter+1; 
                [IF,K] = GlobalK_fast_ver(U,Node,truss,angles);
                R = lmd*F-IF;   MRS = [F,R];
                MUL(FreeDofs,:) = K(FreeDofs,FreeDofs)\MRS(FreeDofs,:);
                dUp = MUL(:,1); dUr = MUL(:,2);
                if iter==1, dUr = 0*dUr; end
                dlmd=nlsmgd(icrm,iter,dUp,dUr,b_lambda);
                dUt = dlmd*dUp+dUr;
                U = U+dUt;
                err = norm(dUt(FreeDofs));
                lmd = lmd+dlmd;
                fprintf('    iter = %d, err = %6.4f, dlambda = %6.4f\n',iter,err,dlmd);
                if err > 1e8, disp('Divergence!'); break; end
            end
    
            if iter>15
                b_lambda = b_lambda/2;
                disp('Reduce constraint radius...')
                icrm = icrm-1;
                U = Uhis(:,max(icrm,1));  % restore displacement
                lmd = Fhis(max(icrm,1));   % restore load
            elseif iter<3
                disp('Increase constraint radius...')
                b_lambda = b_lambda*1.5;
                Uhis(:,icrm) = U;
                Fhis(icrm) = lmd; 
            else
                Uhis(:,icrm) = U;
                Fhis(icrm) = lmd; 
            end
        end
    
    elseif strcmpi(AnalyInputOpt.LoadType, 'Displacement')
        Uhis = zeros(3*size(Node,1),AnalyInputOpt.DispStep*2);
        if isfield(AnalyInputOpt,'Load')
            Fdsp = AnalyInputOpt.Load/AnalyInputOpt.DispStep;
        else
            Fdsp = AnalyInputOpt.AdaptiveLoad(Node,U,1);
        end
        ImpDofs = find(Fdsp~=0);
        FreeDofs = setdiff(setdiff(1:3*size(Node,1),truss.FixedDofs),ImpDofs);
        icrm = 0;  dspmvd = 0;  attmpts = 0;
        mvstepsize = 1;  damping = 1;
        Fhis = zeros(AnalyInputOpt.DispStep,numel(ImpDofs)); 
        while (dspmvd <= 1 && ~AnalyInputOpt.StopCriterion(Node,U,icrm)) && attmpts <= 20
            icrm = icrm+1;
            iter = 0; err = 1;   
            fprintf('icrm = %d, dspimps = %6.4f\n',icrm,dspmvd);
            if isfield(AnalyInputOpt,'AdaptiveLoad')
                Fdsp = AnalyInputOpt.AdaptiveLoad(Node,U,icrm);
            end
            U = U+mvstepsize*Fdsp;
            U(truss.FixedDofs)=0;
            while err>tol && iter<MaxIter
                iter = iter+1;
                [IF,K] = GlobalK_fast_ver(U,Node,truss,angles);
                dU = zeros(3*size(Node,1),1);
                dU(FreeDofs) = K(FreeDofs,FreeDofs)\(-IF(FreeDofs));
                err = norm(dU(FreeDofs));
                U = U+damping*dU; 
                fprintf('    iter = %d, err = %6.4f\n',iter,err);
            end
    
            if iter>=((mvstepsize>1)+1)*MaxIter/(damping+1)  
                % an aggressive step needs more iterations
                attmpts = attmpts+1;
                icrm = icrm-1;
                if attmpts<=10  
                    mvstepsize = mvstepsize*0.5; 
                    disp('Take a more conservative step...')
                else
                    mvstepsize = max(mvstepsize,1)*1.5;  
                    damping = damping*0.75;
                    disp('Take a more aggressive step...')
                end
                U = Uhis(:,max(icrm,1)); % restore displacement            
            else
                dspmvd = dspmvd+mvstepsize/AnalyInputOpt.DispStep;
                attmpts = 0;
                damping = 1;
                if mvstepsize<1
                    mvstepsize = min(mvstepsize*1.1,1); % gradually go back to 1
                else
                    mvstepsize = max(mvstepsize*0.9,1);
                end
                Uhis(:,icrm) = U;
                [Fend,~] = GlobalK_fast_ver(U,Node,truss,angles);
                Fhis(icrm,:) = -Fend(ImpDofs)'; 
            end
        end
    else
        disp('Unknown load type!!!')
    end
    
    icrm = icrm+1;
    Uhis(:,icrm:end) = [];
    Fhis(icrm:end,:) = [];
    end
    
    %--------------------------------------------------------------------------
    function dl=nlsmgd(step,ite,dup,dur,cmp)
    % Modified Generalized Displacement Control Method
    global dupp1 sinal dupc1 numgsp
    if ite==1
        if step==1
            sinal=sign(dot(dup,dup));
            dl=cmp;
            numgsp=dot(dup,dup);   
        else
            sinal=sinal*sign(dot(dupp1,dup));
            gsp=numgsp/dot(dup,dup);
            dl=sinal*cmp*sqrt(gsp);
        end 
        dupp1=dup;
        dupc1=dup;
    else
        dl=-dot(dupc1,dur)/dot(dupc1,dup);
    end
    end