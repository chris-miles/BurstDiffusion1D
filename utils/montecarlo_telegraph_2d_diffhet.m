function particles_at_end = montecarlo_telegraph_2d_diffhet(Tfinal,bdry, Xsource, kon, koff, Dfunc, ron,roff)


dt=1e-3;


[birth_times,~] = generate_birthtimes_telegraph(Tfinal,kon,ron,roff);
birth_times(birth_times>Tfinal) =[];
birth_times = birth_times';

tot_particles_birth = length(birth_times);



survived = ones(1,tot_particles_birth); % vector to store whether they will survive at the end

life_times = exprnd(1/koff,1,tot_particles_birth); % generate lifetimes

death_times = birth_times + life_times; % death is birth + lifetime

survived(death_times<Tfinal) = 0; % if the death time is before T, they die

final_positions = zeros(tot_particles_birth,2); % now we must move the remaining ones to see if they exit nucleus


for n=1:tot_particles_birth
    if survived(n) % if they survived
        steps_to_move = ceil((Tfinal-birth_times(n))/dt); % number of steps to generate

        randn_vals = randn(steps_to_move,2);

        position_n = zeros(steps_to_move,2);
        position_n(1,:) = Xsource;

        for s = 1:steps_to_move
            curr_pos = position_n(s,:);
            new_pos = curr_pos + sqrt(2*dt*Dfunc(curr_pos(1),curr_pos(2)))*randn_vals(s,:);
            position_n(s+1,:) = new_pos;
        end


        exits = find(~insidepoly(position_n(:,1),position_n(:,2),bdry(1,:),bdry(2,:)),1);

        
        %position_n = cumsum(moves_n,1)+Xsource; % add up total moves to the end



        %pos_norm = min([position_n, R-position_n]);


        %sb=sqrt(2*D*dt);
        %escape_probs = exp(-2*(R-pos_norm(1:end-1)).*(R-pos_norm(2:end))/(sb^2));
        %randvals = rand(size(escape_probs));
        %andrewsbray_escape_fix = find(randvals<escape_probs,1);


       % exits = find(~insidepoly(position_n(:,1),position_n(:,2),bdry(1,:),bdry(2,:)),1);
        % if isempty(exits)&&~isempty(andrewsbray_escape_fix)
        %     disp('fix helped');
        %  end

        %if ~isempty(andrewsbray_escape_fix)&&isempty(exits)
        %    disp('andrewsbrayfix mattered')
        %end

        if ~isempty(exits)%||~isempty(andrewsbray_escape_fix) % if it leaves, kill it
            survived(n)=0;
        else % otherwise store its final position
            final_positions(n,:) = position_n(end,:);
            %end
        end
    end
end

% output remaining particles

particles_at_end = final_positions(logical(survived),:);

