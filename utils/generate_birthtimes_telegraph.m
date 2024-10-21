function [birth_times, y] = generate_birthtimes_telegraph(Tfinal,kon,ron,roff)


birth_times = [];
x = rand < ron/(ron+roff);
t = 0;
y=0;


while t <= Tfinal
    r_off = roff*(x==1)*1.0;
    r_on = ron*(x==0)*1.0;
    r_birth = kon*(x==1)*1.0;

    r_tot = r_off+r_on+r_birth;

    % Calculate time to next event
    dt = exprnd(1/r_tot);

    % Update time
    t = t + dt;

    % Determine which event occurs
    event = rand() * r_tot;
    if event < r_off
        x = 0;  %  off
    elseif event < r_off + r_on
        x = 1;  % on
    else 
        birth_times = [birth_times;t];  % birth
        y = y+1;
    end 
end

end
