%
% function Fsdifclose(file)
%
% close sdif file 
%
% INPUT :
%  file  = file handle previously opened with Fsdifopen
%
% OUTPUT :
%  NONE
%
% SEE also : Fsdifopen, Fsdifclose, Fsdifread, and the low level handlers
%     Fsdif_read_handler and Fsdif_write_handler
%
%
% AUTHOR : Axel Roebel
% DATE   : 21.01.2008
%
% $Revision: 1.1 $    last changed $Date: 2008/01/22 00:52:56 $
%
%                                                 Copyright (c) 2008 by  IRCAM 

function Fsdifclose(file)
  if( ~Fsdif_read_handler('close',file) )
    if (~Fsdif_write_handler('close',file) )
      error('warning file handle was not open')
    end
  end
  
  