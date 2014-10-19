%
% Called by setCorePaths
%
% Out:
% - the path to the core in this machine
%
% History:
% 14 March 05 - fpereira - created
%

function [corePath] = getCorePath()

  global HOSTNAME
  % depending on where we are running whomever calls this, the data
  % location will be different
  if isempty(HOSTNAME)
    % can set this globally (recommended if running with jvm)
    % or ask operating system (works fine when 'matlab -nojvm')
    [s,HOSTNAME] = unix('hostname');
    HOSTNAME = lower(deblank(HOSTNAME));
  end;
  
  % check if it is a BIRC machine
  [P,R] = strtok(HOSTNAME,'.');
  if isequal(R,'.c1.imaginginst.edu'); HOSTNAME = 'BIRC'; end
  
  switch HOSTNAME
   case {'human','human.learning.cs.cmu.edu'}
    % Tom's laptop
    corePath = '';
   case {'humanWindows','humanWindows.learning.cs.cmu.edu'}
    % Tom's laptop in Windows2000
    corePath = '';
   case {'pingo','pingo.pc.cs.cmu.edu'}
    % francisco's laptop
    corePath = '/home/fpereira/projecto/FMRI/fmri_core_new';
   case {'gs2011','gs2011.sp.cs.cmu.edu'}
    % rebecca's machine
    corePath = '';
   case {'BIRC'}
    % BIRC machines
    corePath = '/home/fpereira/fmri_core_new';
    
   otherwise
    
    % EVERY MACHINE ON SFS
    corePath = '/afs/cs/project/theo-72/fmri_core_new';
  end
