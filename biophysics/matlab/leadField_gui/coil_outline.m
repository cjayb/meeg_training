function coilframe = coil_outline(vLoc, RotMatrix, coilname)
  %
  % function NMchip = crNMchip(sensLoc, sensNames, units)
  %
  % -----    -----  This is a gradiometer coil.
  % |    |  |    |
  % |a)   \/     | 
  % |     /\     |
  % |    |  | b) |
  % ----- c) -----

  units = 'm';
  
  if strcmp(units, 'mm'),
      unitfactor = 1;
  elseif strcmp(units, 'm'),
      unitfactor = 0.001;
  else
      unitfactor = 0.001;
  end
    
  chipside = 28*unitfactor; %physical (square) chip size in mm
  magside = 26.9*unitfactor;
  glongouts = 24.8*unitfactor;
  gshorts = 9.5*unitfactor;
  gwingsep = 3.2*unitfactor;
  glongins = 10.6*unitfactor;
  
  % chip centre point same for all coils
%   vLoc = SensLoc(1:3);
  % chip rotation matrix
%   RotMatrix = reshape(SensLoc(4:12), 3,3);
  
  % Frame of right gradiometer wing, left its negative
  posframe = [0 0 0;
      gwingsep/2 glongouts/2-glongins 0;
      gwingsep/2 glongouts/2 0;             %2
      gwingsep/2+gshorts glongouts/2 0;     %3
      gwingsep/2+gshorts 0 0;
      gwingsep/2+gshorts -glongouts/2 0;
      gwingsep/2 -glongouts/2 0;
      gwingsep/2 -glongouts/2+glongins 0;
%       gwingsep/2+gshorts/2 glongouts/4 0;
%       gwingsep/2+gshorts/2 0 0;
%       gwingsep/2+gshorts/2 -glongouts/4 0;
      0 0 0];
  
  negframe = -posframe;


  % Combine frames and triangulations of two wings
  frame = [posframe; negframe(2:(length(negframe)),:)];
  
  coilStruct.xgframe = frame;
  % y-gradiometer x-axis is a turn 90 degrees from x-grad x-axis
  Rot = [0 -1 0; 1 0 0; 0 0 1]; % rotation matrix for right-hand rule!!!!
  coilStruct.ygframe = (Rot*frame')';
  
  % Form and triangulate magnetometer
  coilStruct.magframe = [magside/2 magside/2 0;...
		    magside/2 -magside/2 0;...
		    -magside/2 -magside/2 0;...
		    -magside/2 magside/2 0;...
		    magside/2 magside/2 0];
  		     
  % Form chip outline for plotting
  coilStruct.outline = [chipside/2 chipside/2 0;...
		    chipside/2 -chipside/2 0;...
		    -chipside/2 -chipside/2 0;...
		    -chipside/2 chipside/2 0;...
		    chipside/2 chipside/2 0];
  
  
  coilStruct.xgframe = (RotMatrix*coilStruct.xgframe')' ...
      + repmat(vLoc, size(coilStruct.xgframe, 1), 1);
  coilStruct.ygframe = (RotMatrix*coilStruct.ygframe')' ...
      + repmat(vLoc, size(coilStruct.ygframe, 1), 1);
  coilStruct.magframe = (RotMatrix*coilStruct.magframe')' ...
      + repmat(vLoc, size(coilStruct.magframe, 1), 1);
  coilStruct.outline = (RotMatrix*coilStruct.outline')' ...
      + repmat(vLoc, size(coilStruct.outline, 1), 1);

  

  % save location of chip center point
  coilStruct.loc = vLoc;
  
  % save rotation matrix
  coilStruct.RotMatrix = RotMatrix;

  if strcmp(coilname, 'mag')
    coilframe = coilStruct.magframe;
  elseif strcmp(coilname, 'xgrad')
    coilframe = coilStruct.xgframe;
  elseif strcmp(coilname, 'ygrad')
    coilframe = coilStruct.ygframe;
  elseif strcmp(coilname, 'outline')
    coilframe = coilStruct.outline;
  end
    