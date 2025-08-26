% function rgb = my_color(color_name)
% 
% color_name:
% 
% 'light_green'
% 'maroon'    
% 'pale_blue' 
% 'butterscotch'       
% 'indigo'         
% 'pink'       
% 'yellow'       
% 'orange'        
% 'green'         
% 'blue'      
% 'magenta' 
% 'light_orange'  
% 'pale_orange' 
% 'bright_green'   
% 'bright_blue'
% 'light_red'
% 'light_blue'
% 'grey'
% 'dark_green'

function rgb = my_color(color_name)

rgb = [0 0 0];

if strcmp(color_name, 'light_green')
    rgb = [51,204,0];
elseif strcmp(color_name, 'maroon')
    rgb = [204,0,62];
elseif strcmp(color_name, 'pale_blue')
    rgb = [92,181,229];
elseif strcmp(color_name, 'butterscotch')
    rgb = [246,167,53];
elseif strcmp(color_name, 'indigo')
    rgb = [84,12,106];
elseif strcmp(color_name, 'pink')
    rgb = [249,3,203];
elseif strcmp(color_name, 'yellow')
    rgb = [233,255,36];
elseif strcmp(color_name, 'orange')
    rgb = [235,91,19];
elseif strcmp(color_name, 'green')
    rgb = [0,150,9];
elseif strcmp(color_name, 'blue')
    rgb = [91,19,235];
elseif strcmp(color_name, 'magenta')
    rgb = [237,31,179];
elseif strcmp(color_name, 'light_orange')
    rgb = [255,197,7];
elseif strcmp(color_name, 'pale_orange')
    rgb = [240,117,78];
elseif strcmp(color_name, 'bright_green')
    rgb = [35,231,52];
elseif strcmp(color_name, 'bright_blue')
    rgb = [7,205,255];
elseif strcmp(color_name, 'light_red')
    rgb = [236,214,214];
elseif strcmp(color_name, 'light_blue')
    rgb = [205,227,247];
elseif strcmp(color_name, 'grey')
    rgb = [.8,.8,.8]*256;
elseif strcmp(color_name, 'dark_green')
    rgb = [0,.5,0]*256;
end
rgb = rgb/256;    