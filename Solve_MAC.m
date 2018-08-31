if ~isfield(display,'calculated') || ~isfield(display.calculated,'figure')
    display.calculated.figure = figure();
    if display.set_Figure_Position &&~ video.record
        display.calculated.figure.Position=display.figure_Position;
    end
else
    if ~isgraphics(display.calculated.figure,'figure')
        display.calculated.figure = figure();
        if display.set_Figure_Position &&~ video.record
            display.calculated.figure.Position=display.figure_Position;
        end
    end
end
lastwarn('')%clearing warnings...


time_step = 0;
pde.calculated.u0 = @(x,y) pde.u0t(x,y,0);
pde.calculated.v0 = @(x,y) pde.v0t(x,y,0);
pde.calculated.un = @(x,y) pde.unt(x,y,0);
pde.calculated.vn = @(x,y) pde.vnt(x,y,0);
pde.calculated.fu = @(x,y) pde.fut(x,y,0);
pde.calculated.fv = @(x,y) pde.fvt(x,y,0);


if pde.dynamic.dynamic && video.record
    if display.full_HD_for_video
       factor = 1+display.monitor_4k;
       display.calculated.figure.Position = [(factor-1)*1920/4 (factor-1)*1088/4 1920/factor 1088/factor];
    end

    try
        VideoObj = VideoWriter(video.name,video.profile); 
        VideoObj.FrameRate = video.FrameRate; 
        open(VideoObj); 
    catch  ME
        switch ME.identifier
            case 'MATLAB:audiovideo:VideoWriter:fileNotWritable'
                warning('Videofile was not writable. Trying to close it...')
                close(VideoObj);
                VideoObj = VideoWriter(video.name,video.profile); 
                VideoObj.FrameRate = video.FrameRate;
                open(VideoObj);
                disp('Successfully reopend the Videofile');
            otherwise
                rethrow(ME)
        end
    end
end

domain.calculated.on_domain = @(x,y) x>=0 & y>=0 & x<=domain.lx & y<=domain.ly &~domain.on_obstacle(x,y);

display.calculated.hx = domain.lx/display.nx;
display.calculated.hy = domain.ly/display.ny;


mesh.calculated = Generate_Mesh(domain,mesh,display);

                                                        
%Compute System Matrix M:
UxU = Create_UxU_Block(mesh);
%figure();spy(UxU);

UxV = sparse(mesh.calculated.nu,mesh.calculated.nv);

UxP = Create_UxP_Block(mesh);


VxV = Create_VxV_Block(mesh);
%figure();spy(VxV);

VxP = Create_VxP_Block(mesh);

PxP = sparse(mesh.calculated.np,mesh.calculated.np);

M = [  UxU , UxV , UxP;
       UxV', VxV , VxP;
       UxP', VxP', PxP];

if(debugging.cleanup)
   clear UxU UxV UxP VxV VxP PxP
end
   
%figure();spy(M);

%Apply Forces
F = sparse(mesh.calculated.DOFs,1);
F(1:mesh.calculated.nu) = pde.calculated.fu(mesh.calculated.X_u(:)',mesh.calculated.Y_u(:)');
F(mesh.calculated.nu+1:mesh.calculated.nu+mesh.calculated.nv) = pde.calculated.fv(mesh.calculated.X_v(:)',mesh.calculated.Y_v(:)');


[M,F] = Restrict_To_Domain(M,F,mesh.calculated.domain_vec);


[M,F] = Apply_BCs_U( M,F,domain,mesh,pde);
[M,F] = Apply_BCs_V( M,F,domain,mesh,pde);

for i=find(sum(M(1:mesh.calculated.nu+mesh.calculated.nv,mesh.calculated.nu+mesh.calculated.nv+(1:mesh.calculated.np))')~=0) 
   %find missing neighboring p's
   M(i,mesh.calculated.nu+mesh.calculated.nv+(1:mesh.calculated.np))=0; 
end


meanline = mesh.h^2*[zeros(1,mesh.calculated.nu+mesh.calculated.nv),ones(1,mesh.calculated.np)]; %a matrix line representing the pressure mean


if debugging.rankcheck
    Mfull = full(M);

    disp(['Rank of Mfull = ', num2str(rank(Mfull)), ' out of ', num2str(np+nu+nv),' DOFs'])
    disp(['Rank of u submatrix = ', num2str(rank(Mfull(1:nu,1:nu))),' out of ', num2str(nu)])
    disp(['Rank of v submatrix = ', num2str(rank(Mfull(nu+(1:nv),nu+(1:nv)))),' out of ', num2str(nv)])

    if rank(Mfull) <nu+nv+np-1
       return 
    end
end

if ~pde.dynamic.dynamic

    u_sol = [M;meanline]\[F;0];
    %u_sol = lsqminnorm([M;meanline],[F;0]);

    % reordering
    U = reshape(u_sol(1:mesh.calculated.nu),mesh.calculated.nx+1,mesh.calculated.ny)';
    V = reshape(u_sol(mesh.calculated.nu+1:mesh.calculated.nu+mesh.calculated.nv),mesh.calculated.nx,mesh.calculated.ny+1)';
    P = reshape(u_sol(mesh.calculated.nu+mesh.calculated.nv+1:end),mesh.calculated.nx,mesh.calculated.ny)';


    %% Plotting routine
    Plot_Solution; %TODO FUNCTION!
else
    u_sol = zeros(mesh.calculated.DOFs,1);
    u_start_matrix = pde.dynamic.u_start(mesh.calculated.X_u,mesh.calculated.Y_u);
    u_sol(1:mesh.calculated.nu)=u_start_matrix(:);
    v_start_matrix = pde.dynamic.v_start(mesh.calculated.X_v,mesh.calculated.Y_v);
    u_sol(mesh.calculated.nu+(1:mesh.calculated.nv))=v_start_matrix(:);
    u_sol(~mesh.calculated.domain_vec) = 0;
    
    S = spdiags([ones(mesh.calculated.nu+mesh.calculated.nv,1);zeros(mesh.calculated.np,1)],0,mesh.calculated.DOFs,mesh.calculated.DOFs);
    S(mesh.calculated.on_bndry_u,:) = 0;
    S(mesh.calculated.on_bndry_v+mesh.calculated.nu,:)=0;
    %u_stat = M\F;%stationary solution
    %u_sol(nu+nv+(1:np))=u_stat(nu+nv+(1:np));
    
    %P_Mat = M(1:nu+nv,nu+nv+(1:np));
    
    
    %p = P_Mat\(F(1:nu+nv)-M(1:nu+nv,1:nu+nv)*u_sol(1:nu+nv));
    %u_sol(nu+nv+(1:np))=p;
    
    
    
    %u_sol(on_bndry_u) = u0(X_u(on_bndry_u),Y_u(on_bndry_u));
    %u_sol(on_bndry_v+nu) =v0(X_v(on_bndry_v),Y_v(on_bndry_v));
    
    integrator.calculated.dT = integrator.T/integrator.nT;
    
    
    
    U = reshape(u_sol(1:mesh.calculated.nu),mesh.calculated.nx+1,mesh.calculated.ny)';
    V = reshape(u_sol(mesh.calculated.nu+1:mesh.calculated.nu+mesh.calculated.nv),mesh.calculated.nx,mesh.calculated.ny+1)';
    P = reshape(u_sol(mesh.calculated.nu+mesh.calculated.nv+1:end),mesh.calculated.nx,mesh.calculated.ny)';
    Plot_Solution;
    if video.record
        frame = getframe(display.calculated.figure);
        writeVideo(VideoObj, frame);
    end
    for time_step=1:integrator.nT
        %u_sol = (dT*theta*M+speye(size(M)))\(u_sol + dT*theta*F+dT*(1-theta)*(F-M*u_sol));
        %u_sol = u_sol+dT*(F-M*u_sol);%explizit euler
        %u_sol = (dT*M+speye(size(M)))\(F+u_sol);%implizit euler
        
        
        if pde.dynamic.time_Dependent
            pde.calculated.fu = @(x,y) pde.fut(x,y,time_step*integrator.calculated.dT);
            pde.calculated.fv = @(x,y) pde.fvt(x,y,time_step*integrator.calculated.dT);
            F(1:mesh.calculated.nu) = pde.calculated.fu(mesh.calculated.X_u(:)',mesh.calculated.Y_u(:)');
            F(mesh.calculated.nu+1:mesh.calculated.nu+mesh.calculated.nv) = pde.calculated.fv(mesh.calculated.X_v(:)',mesh.calculated.Y_v(:)');
            pde.calculated.u0 = @(x,y) pde.u0t(x,y,time_step*integrator.calculated.dT);
            pde.calculated.v0 = @(x,y) pde.v0t(x,y,time_step*integrator.calculated.dT);
            [M,F] = Apply_BCs_U( M,F,domain,mesh,pde);
            [M,F] = Apply_BCs_V( M,F,domain,mesh,pde);
        end
        
        
             
        u_sol = [(S+integrator.calculated.dT*integrator.theta*M);integrator.calculated.dT*integrator.theta*meanline]\...
            [((S+integrator.calculated.dT*(1-integrator.theta)*M)*u_sol+integrator.calculated.dT*F);0];
        if debugging.stop_video_on_rank_deficit
            [~,warn] = lastwarn;
            if strcmp(warn,'MATLAB:nearlySingularMatrix') || strcmp(warn,'MATLAB:rankDeficientMatrix')
                warning('The time evolution Matrix is singular. In order to resume anyway, please set stop_video_on_rank_deficit to false')
                break;
            end
        end
        
        %u_sol(on_bndry_u) = u0(X_u(on_bndry_u),Y_u(on_bndry_u));
        %u_sol(on_bndry_v+nu) =v0(X_v(on_bndry_v),Y_v(on_bndry_v));
    
        
        %PLOTTING!!!
        if ~mod(time_step,integrator.plot_every)
            U = reshape(u_sol(1:mesh.calculated.nu),mesh.calculated.nx+1,mesh.calculated.ny)';
            V = reshape(u_sol(mesh.calculated.nu+1:mesh.calculated.nu+mesh.calculated.nv),mesh.calculated.nx,mesh.calculated.ny+1)';
            P = reshape(u_sol(mesh.calculated.nu+mesh.calculated.nv+1:end),mesh.calculated.nx,mesh.calculated.ny)';
            Plot_Solution;
            drawnow();
            if video.record
                frame = getframe(display.calculated.figure);
                writeVideo(VideoObj, frame);
            end
        end
    end
    
end

if pde.dynamic.dynamic && video.record
    close(VideoObj); 
end

