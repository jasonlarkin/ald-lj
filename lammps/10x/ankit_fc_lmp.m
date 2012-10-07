
%--------------------------------------------------------------------------
%-------------------------INPUT--------------------------------------------
%--------------------------------------------------------------------------
LD.constant.kb = 1.3806E-23;                    %aJ/k (1.3806E-23 J/K)
LD.constant.hbar = 1.054E-34;                %J/s
LD.constant.i = sqrt(-1);
LD.constant.c = 29979245800.00019;      %cm/s
LD.constant.s2ps = 10E-13;
LD.constant.ang2m = 10E-11;
LD.constant.eV2J = 1.60217E-19;

format long

%--------------------------------------------------------------------------
%-------------------------READ---------------------------------------------
%--------------------------------------------------------------------------

LD.FC2.id = load('./FULLorder2.dat');
% LD.FC3.id = load('./update_list/FULLorder3_2.dat');

LD.pos = load('./FULLpos.dat');

%--------------------------------------------------------------------------
%LD Input
%--------------------------------------------------------------------------
%define supercell integers
LD.Nx = 6;                LD.Ny = 6;                LD.Nz = 6;

%cell parameters
LD.cell = [1.56 1.56 1.56 90 90 90 0 0 0 0 0 0];
%lattice vectors
LD.latvec =...
        [0.0 0.5 0.5;
		0.5 0.0 0.5;
		0.5 0.5 0.0];
LD.latvec = LD.latvec*LD.cell(1);
LD.latvec_rec = [-1.0 1.0 1.0;
			1.0 -1.0 1.0;
			1.0 1.0 -1.0];
        
% LD.x.ucell.frac =...
%     [0.1250000000000000  0.1250000000000000  0.1250000000000000
%     0.8750000000000000  0.8750000000000000  0.8750000000000000];
%                 
% LD.x.ucell.gulp =...
%     [0.8750000000000000  0.8750000000000000  0.8750000000000000 0 1 1 1
%     0.1250000000000000  0.1250000000000000  0.1250000000000000 0 1 1 1];


plot3(LD.pos(:,1),LD.pos(:,2),LD.pos(:,3),'.')

%--------------------------------------------------------------------------
%LD Input
%--------------------------------------------------------------------------

LD.pos(:,2:4) = LD.pos(:,2:4)*LD.cell(1); 

%--------------------------------------------------------------------------
%Prepare LAMMPS Data
%--------------------------------------------------------------------------
lammps.alpha = (180/pi)* ...
    atan2(norm(cross(LD.latvec(1,:),LD.latvec(2,:))),...
    dot(LD.latvec(1,:),LD.latvec(2,:)));
lammps.beta = lammps.alpha;
lammps.gamma = lammps.beta;

lammps.lx = LD.cell(1)*LD.Nx;
lammps.xy = LD.cell(2)*cos(LD.cell(6));
lammps.xz = LD.cell(3)*cos(LD.cell(5));
lammps.ly = sqrt(LD.cell(2)^2 - lammps.xy^2);
lammps.yz = ( LD.cell(2)*LD.cell(3)*cos(LD.cell(4)) - ...
    lammps.xy*lammps.xz ) / lammps.ly;
lammps.lz = sqrt( LD.cell(3)^2 - lammps.xz^2 - lammps.yz^2 );

%--------------------------------------------------------------------------
%pause
%--------------------------------------------------------------------------

%set the lammps ids
    LD.id(1:size(LD.pos,1),1) = 1:size(LD.pos,1); 
    LD.type(1:size(LD.pos,1),1) = 1;
%initialize  
    LD.dx = 0.0001;
    LD.precision = int2str(log10(ceil(1/LD.dx))+4);
    
    LD.pos_tmp = LD.pos(:,:);    
    LAMMPS.force0 = ankit_lmp_force(LD);
    
%--------------------------------------------------------------------------
%pause  
%-------------------------------------------------------------------------- 

%--------------------------------------------------------------------------
tic
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%FC2
%--------------------------------------------------------------------------

for iFC = 1:1:size(LD.FC2.id,1)
%initialize    
    LD.pos_tmp = LD.pos(:,:);
%increment position for iFC,2
    LD.pos_tmp( LD.FC2.id(iFC,1) +1 , LD.FC2.id(iFC,3) + 1 ) = ...
        LD.pos_tmp( LD.FC2.id(iFC,1) +1 , LD.FC2.id(iFC,3) + 1 ) + LD.dx;
%calc changed force    
    LAMMPS.force = ankit_lmp_force(LD);
%change in force on iFC,1
    LD.FC2.phi(iFC) =...
        -(LAMMPS.force(LD.FC2.id(iFC,2) +1 , LD.FC2.id(iFC,4)  ) - ...
        LAMMPS.force0(LD.FC2.id(iFC,2) +1 , LD.FC2.id(iFC,4) ) ) / LD.dx ;
%--------------------------------------------------------------------------
%pause  
%--------------------------------------------------------------------------  
end

plot(LD.FC2.phi)
    
%--------------------------------------------------------------------------
%pause  
%--------------------------------------------------------------------------     
    
%--------------------------------------------------------------------------
%FC3
%--------------------------------------------------------------------------

% for iFC = 1:1:size(LD.FC3.id,1)
%     
% %initialize    
%     LD.pos_tmp = LD.pos(:,:);
% %increment position beta
%     LD.pos_tmp(LD.FC3.id(iFC,2) , LD.FC3.id(iFC,5) ) = ...
%         LD.pos_tmp(LD.FC3.id(iFC,2) , LD.FC3.id(iFC,5) ) - LD.dx;
% %increment position gamma
%     LD.pos_tmp(LD.FC3.id(iFC,3) , LD.FC3.id(iFC,6) ) = ...
%         LD.pos_tmp(LD.FC3.id(iFC,3) , LD.FC3.id(iFC,6) ) - LD.dx;
% %calc changed force
%     LAMMPS.force1 = ankit_lmp_force(LD);
% 
% %initialize    
%     LD.pos_tmp = LD.pos(:,:);
% %increment position beta
%     LD.pos_tmp(LD.FC3.id(iFC,2) , LD.FC3.id(iFC,5) ) = ...
%         LD.pos_tmp(LD.FC3.id(iFC,2) , LD.FC3.id(iFC,5) ) - LD.dx;
% %increment position gamma
%     LD.pos_tmp(LD.FC3.id(iFC,3) , LD.FC3.id(iFC,6) ) = ...
%         LD.pos_tmp(LD.FC3.id(iFC,3) , LD.FC3.id(iFC,6) ) + LD.dx;
% %calc changed force
%     LAMMPS.force2 = ankit_lmp_force(LD);
%     
% %initialize    
%     LD.pos_tmp = LD.pos(:,:);
% %increment position beta
%     LD.pos_tmp(LD.FC3.id(iFC,2) , LD.FC3.id(iFC,5) ) = ...
%         LD.pos_tmp(LD.FC3.id(iFC,2) , LD.FC3.id(iFC,5) ) + LD.dx;
% %increment position gamma
%     LD.pos_tmp(LD.FC3.id(iFC,3) , LD.FC3.id(iFC,6) ) = ...
%         LD.pos_tmp(LD.FC3.id(iFC,3) , LD.FC3.id(iFC,6) ) - LD.dx;
% %calc changed force
%     LAMMPS.force3 = ankit_lmp_force(LD);
%     
% %initialize    
%     LD.pos_tmp = LD.pos(:,:);
% %increment position beta
%     LD.pos_tmp(LD.FC3.id(iFC,2) , LD.FC3.id(iFC,5) ) = ...
%         LD.pos_tmp(LD.FC3.id(iFC,2) , LD.FC3.id(iFC,5) ) + LD.dx;
% %increment position gamma
%     LD.pos_tmp(LD.FC3.id(iFC,3) , LD.FC3.id(iFC,6) ) = ...
%         LD.pos_tmp(LD.FC3.id(iFC,3) , LD.FC3.id(iFC,6) ) + LD.dx;
% %calc changed force
%     LAMMPS.force4 = ankit_lmp_force(LD);
% 
%     
% %change in force on alpha
%     LD.FC3.phi(iFC) =...
%         -...
%         (...
%         LAMMPS.force4(LD.FC3.id(iFC,1) , LD.FC3.id(iFC,4) ) - ...
%         LAMMPS.force3(LD.FC3.id(iFC,1) , LD.FC3.id(iFC,4) ) - ...
%         LAMMPS.force2(LD.FC3.id(iFC,1) , LD.FC3.id(iFC,4) ) + ...
%         LAMMPS.force1(LD.FC3.id(iFC,1) , LD.FC3.id(iFC,4) ) ...
%         )...
%         / (4*LD.dx^2);
% %--------------------------------------------------------------------------
% %pause  
% %-------------------------------------------------------------------------- 
% 
% end
% 
% plot(LD.FC3.phi)

%--------------------------------------------------------------------------
%pause  
%-------------------------------------------------------------------------- 


%--------------------------------------------------------------------------
%pause  
%-------------------------------------------------------------------------- 

lj = m_lj;

LD.FC2.phi = LD.FC2.phi*lj.eps/lj.sigma / LD.constant.eV2J * LD.constant.ang2m;


dlmwrite('./PHI2.dat',LD.FC2.phi','delimiter',' ','precision', '%10.5f');
% dlmwrite('./PHI3.dat',LD.FC3.phi','delimiter',' ','precision', '%10.5f');



