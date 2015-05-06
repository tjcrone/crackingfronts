function [] = progressbar(n,totn,mfilename,textinfo)
%this generates a command line progress bar. n must increment by 1 and
%start with 1

olddone = floor(100*(n-1)/totn);
done = floor(100*n/totn);

if olddone == done
else

    %compute remaining time
    timespent = toc;
    timeperstep = timespent/n;
    remtime = round((totn-n)*timeperstep);
    h = floor(remtime/60/60);
    m = floor(remtime/60)-60*h;
    s = remtime-60*60*h-60*m;

    %setup strout
    done50 = floor(n/totn*100/2);
    left50 = 50-done50;
    strout1 = ['[' repmat('=',1,done50) repmat(' ',1,left50) ']'];
    %strout2 = [strout1 sprintf(' %3.0f%% %3.0fh%2.0fm%2.0fs',done,h,m,s) textinfo];
    %stroutlength = length(strout2);

    %clear line
    %if n~=1
    %    fprintf(1,repmat('\b',1,stroutlength));
    %end

    %clear screen
    clc;

    %print strout
    fprintf(1,'%s %s\n',mfilename,textinfo);
    fprintf(1,strout1);
    fprintf(1,' %3.0f%% %3.0fh%2.0fm%2.0fs',done,h,m,s);

    if n==totn
        fprintf(1,'\n');
    end


end
