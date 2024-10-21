function particles_at_end = montecarlo_telegraph_robin(Tfinal,R, Xsource, kon, koff, D, ron,roff,kappa)

dt=1e-4;


[birth_times,~] = generate_birthtimes_telegraph(Tfinal,kon,ron,roff);
birth_times(birth_times>Tfinal) =[];
birth_times = birth_times';

tot_particles_birth = length(birth_times);


survived = ones(1,tot_particles_birth); % vector to store whether they will survive at the end

life_times = exprnd(1/koff,1,tot_particles_birth); % generate lifetimes

death_times = birth_times + life_times; % death is birth + lifetime

survived(death_times<Tfinal) = 0; % if the death time is before T, they die

final_positions = zeros(tot_particles_birth,1); % now we must move the remaining ones to see if they exit nucleus


for n=1:tot_particles_birth
    if survived(n) % if they survived
        steps_to_move = ceil((Tfinal-birth_times(n))/dt); % number of steps to generate
        moves_n = sqrt(2*D*dt)*randn(steps_to_move,1);

        position_n = Xsource; % add up total moves to the end
        exits=0;

        for s=1:steps_to_move
            new_pos = position_n + moves_n(s);
            if abs(new_pos)>R
                escape_event = rand < kappa*sqrt(dt)*sqrt(pi)/sqrt(D);
                if escape_event
                    exits = 1;
                else
                    if new_pos > R
                        new_pos = 2*R-new_pos;
                    else
                        new_pos = -2*R-new_pos;
                    end
                end
            end
            position_n = new_pos;
        end

        %pos_norm = min([position_n, R-position_n]);

        %sb=sqrt(2*D*dt);
        %escape_probs = exp(-2*(R-pos_norm(1:end-1)).*(R-pos_norm(2:end))/(sb^2));
        %randvals = rand(size(escape_probs));
        %andrewsbray_escape_fix = find(randvals<escape_probs,1);


        %exits = find(abs(position_n)>R,1);  % check if the particle movement is ever NOT in the mask

        % if isempty(exits)&&~isempty(andrewsbray_escape_fix)
        %     disp('fix helped');
        %  end

        %if ~isempty(andrewsbray_escape_fix)&&isempty(exits)
        %    disp('andrewsbrayfix mattered')
        %end

        if exits%||~isempty(andrewsbray_escape_fix) % if it leaves, kill it
            survived(n)=0;
        else % otherwise store its final position
            final_positions(n,:) = position_n;
        end
    end
end

% output remaining particles

particles_at_end = final_positions(logical(survived));

