function sim = DesignSimCIGSRunSim(sim)
% Design the simulation

sim.optical.setup = struct(...
     ... % Optical Parameters
    'optical_toggle', 1, ... % Toggle: 1 to run optical simulation, 0 to skip
    'BLL',0, ... % 1 for Beer--Lambert law and 0 for RCWA 
    'Ez_A_or_iE', 0, ...    % Choice for rebuilding Ez from Hy, 0 for A, 1 for iE
    'pol', 0, ...               % 0: Unpolarized, 1: p-polarized, 2: s-polarized
    ... % Analysis Options
    'RebuildFields', 0, ... % Use ifft to rebuild fields (overridden if Eplot or Hplot toggled)
    ... % Display options
    'optplot', 0, ...         % Produce plots (0: no, 1: yes)
    'epsplot', 0, ... Plot permittivity
    'Eplot', 0, ... Plot Efield
    'Hplot', 0, ... Plot Hfield
    ... % Mesh for angle  (following is untested)
    'degtheta0', 0, ...     % Minimum angle of Incident Light 
    'degtheta1', 0, ...     % Maximum angle of Incident Light 
    'ntheta', 1, ...        % Number of angles 
    ... % Mesh for wavelength
    'nmlambda0', 301, ...   % Minimum wavelength in nm
    'nmlambda1', 1200, ...  % Maximum wavelength in nm
    'nlambda', 1 + floor((1200-301)/10), ...     % Number of wavelengths
    ... % discretization along x direction 
    'Nx', 100, ...       % Number sample points in x-direction (for plotting)
    'nslices', [1, 1,1,100,35,250,1,1, 1], ... % Number of slices in each sec (DesignJunction.m), uses 1 if not set.
    ... % Parameters for Nt convergence 
... % Do we ever use this?  No!!
    'CheckNtConvergence', 0, ... Use adaptive Nt (Use at own risk)
    'Nttol', 1e9,...      % Required tolerance in mA/(nm cm^2) for increasing Nt by 2
    'Nttol_lower', 1e-6, ...%  Required tolerance in mA/(nm cm^2) for decreasing Nt by 2 (used in serial calcs)
    'minNt', 10, ...         % Start value for Nt (tom:40)  %  number of Floquet harmonics for RCWA
    'maxNt', 40, ...        % Largest value for Nt (tom:40)
    'nosuccreq', 2, ... % Number of successes required
    'faillimit', 10, ... % Number of times worse before fail
    'ptol', 0.1 ... % Predictive stepping toler
);
    
sim.electrical.setup = struct(...
    ... % Electrical Parameters
    'electrical_toggle', 1, ... % Toggle: 1 to run electrical simulation, 0 to skip
    'VdV', 0.01, ... % Change in volrage 
    ... % Simulation meshing
    'nslices', [1, 1,1,100,35,250,5,5, 1], ... % Number of slices in each sec (DesignJunction.m
    'nx_sec', [100,35,250], ... % Number of mesh points in each section (electrical)
    ...%Use FD or HDG method
    'FD', 1 , ... %Use the FD method(1) or 0 to use HDG
   ... % Start HDG parameters
    'pdeg',5, ... % Polynomial degree
... % what is this for  (take out)
    'rerunbackward', 0, ... % Recalculate JV backward to check convergence via hystersis
 ... # remove?
    'nmdz', 2, ... % Approx. resolution of mesh (disabled)
    'meshfn' , '1', ...& sin(pi*x)+0.01 , ... % Function to specify relative lengths of each segment, scales to match Lz. 1 for uniform
    ... % Global polys. The poly ax^2+bx+c -> [a, b, c], for example.
    'upwind', 1, ...  % If toggled, will ignore value of t and upwind (HDG only)
    'tensor', 0, ... % Uses tensors internally
    't_1', 1e-6, ...
    't_2', 1e3, ...
    'quadtoggle', 1, ... Toggle quad vs lobatto integration
    'intdeg', 10, ... % Integration accuracy
    ... % Result plots
    'SolverSteps', 0, ... % Display solver steps
    'plotsolve', 0, ... % Number of steps to plot results. 0 for no output. Inf for first and last.
    'ExportSCdetails', 0, ... Export the details (bandgaps, densities, recombination etc.) at short-circuit (untested)
    'nx_plot', 2, ...
    'plot_densities', 0, ...
    'plot_fields', 0, ...
    'plot_levels', 0, ...
    'plot_currents', 0, ...
    'plot_U', 0, ...
    'plot_tau', 0, ...
  ... % Tolerance
    'dampmin', 0.01, ... % Set the initial damping level
    'damprelax', 1.1,... % Rate of relaxation of damping
    'trydamping', 100, ... % >=1 Number of times damping is increases to aid convergence
    'dampingincrease', 5, ... Amount to increase damping by each attempt
    'mesh_increase', 1.5, ... mesh multiplier on fail
    'refine_on_negative', 0, ...
    'refine_on_Jnoise', 1, ...
    'maxJnoiserefinefrac', 0.2, ...
    'refine_always', 0, ...
    'maxloop', 40, ... % Maximum number of iteration for convergence
    'maxpdeg', 12, ... %
    'maxtime', 60*60, ... % Maximum time allowed for the initial HDG convergence
    'mindV', 1e-6, ...
    'reltol', 1e5 ,... % Required relative tolerance
    'abstol', 1e-1, ... % Required absolute tolerance   
    'warn', 1, ... % Display warning messages
    'Jtol', 0.4, ... Standard is 0.1. Larger values allow more noisy J profile.
... % do we do this
    'forcepositive', 0, ... % Force the carriers to be positive
    'test_J', 0, ... %Test the Jacobian
    'PVrefinemax', 100, ... % Max refinement steps allowed
    'PVtol', 1e-3 ... % Tolerance for finding Pmax
);


%%%%%%%%%%%%%%%%%%%%%%
% Do not edit below this box %
%%%%%%%%%%%%%%%%%%%%%%

% Create timer
sim.electrical.setup.start_time = tic;
sim.electrical.setup.total_time = toc(sim.electrical.setup.start_time);

% Make sure nslices is correct lenght
diff =  sim.material.nsec - length(sim.optical.setup.nslices);
if diff > 0
    sim.optical.setup.nslices = padarray(sim.optical.setup.nslices, [0,1], diff, 'post');
end

% Set Nt to minNt
sim.optical.setup.Nt = sim.optical.setup.minNt;

% Total step counter
sim.electrical.setup.counter = 0;
end



