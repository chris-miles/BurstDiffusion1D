function outvol = compute_volume(results)
%https://mail.mathpages.com/home/kmath393/kmath393.htm

elem = results.Mesh.Elements;
sol = results.NodalSolution;
nodes=  results.Mesh.Nodes;

nElem = length(elem);

if isa(results,'pde.StationaryResults')


    %nNodes = length(nodes);

    outvol1  = 0;
    outvol2  = 0;

    for e = 1:nElem
        e1 = elem(1,e);
        e2 = elem(2,e);
        e3 = elem(3,e);

        n1 = nodes(:,e1);
        n2 = nodes(:,e2);
        n3 = nodes(:,e3);

        x1 = n1(1);
        y1 = n1(2);

        x2 = n2(1);
        y2 = n2(2);

        x3 = n3(1);
        y3 = n3(2);

        z1 = sol(e1,1);
        z2 = sol(e2,1);
        z3 = sol(e3,1);



        ve = (z1+z2+z3)*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)/6;


        outvol1 = outvol1+ve;

    end
    outvol = outvol1;

else
    solTimes = results.SolutionTimes;

    nTimes = length(solTimes);
    outvol = zeros(1,nTimes);

    for t=1:nTimes

        solt = sol(:,:,t);

        outvolt1  = 0;

        for e = 1:nElem
            e1 = elem(1,e);
            e2 = elem(2,e);
            e3 = elem(3,e);

            n1 = nodes(:,e1);
            n2 = nodes(:,e2);
            n3 = nodes(:,e3);

            x1 = n1(1);
            y1 = n1(2);

            x2 = n2(1);
            y2 = n2(2);

            x3 = n3(1);
            y3 = n3(2);

            z1 = solt(e1,1);
            z2 = solt(e2,1);
            z3 = solt(e3,1);

            w1 = solt(e1,2);
            w2 = solt(e2,2);
            w3 = solt(e3,2);

            ve = (z1+z2+z3)*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)/6;

            outvolt1 = outvolt1+ve;

        end

        outvol(1,t)=outvolt1;

    end

end

end


