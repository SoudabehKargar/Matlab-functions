
close all

y1 = [1:6];
x1 = 1*ones(1,6);

y2 = [7:9]
x2 = 1*ones(1,3)

y3 = [10:14]
x3 = 1*ones(1,5)

x4 = 2*ones(1,6)
x5 = 2*ones(1,3)
x6 = 2*ones(1,5)

x   = ones(1,14);
xx  = 2*ones(1,14);
y = [1:14];


figure
plot(x1,y1,'ko'), hold on
plot(x2,y2,'k.')
plot(x3,y3,'ko')

plot(x4,y1,'ko'), hold on
plot(x5,y2,'k.')
plot(x6,y3,'ko')

for ii = 1:6
    for kk=1:6
        for jj = 1:1
            plot([jj,jj+1],[kk,ii],'k-'), hold on
        end
% plot([1,3],[1,1],'-')
    end
end

for ii = 10:14
    for kk=10:14
        for jj = 1:1
            plot([jj,jj+1],[kk,ii],'k-'), hold on
        end
% plot([1,3],[1,1],'-')
    end
end

for ii = 1:6
    for kk=10:14
        for jj = 1:1
            plot([jj,jj+1],[kk,ii],'k-'), hold on
        end
% plot([1,3],[1,1],'-')
    end
end

for ii = 10:14
    for kk=1:6
        for jj = 1:1
            plot([jj,jj+1],[kk,ii],'k-'), hold on
        end
% plot([1,3],[1,1],'-')
    end
end

for ii = 1:6
    for kk=6.5:6.5
        for jj = .5:.5
            plot([jj,jj+.5],[kk,ii],'k-'), hold on
        end
    end
end

for ii = 10:14
    for kk=6.5:6.5
        for jj = .5:.5
            plot([jj,jj+.5],[kk,ii],'k-'), hold on
        end
    end
end

for ii = 1:6
    for kk=10.5:10.5
        for jj = .5:.5
            plot([jj,jj+.5],[kk,ii],'k-'), hold on
        end
    end
end

for ii = 10:14
    for kk=10.5:10.5
        for jj = .5:.5
            plot([jj,jj+.5],[kk,ii],'k-'), hold on
        end
    end
end

for ii = 1:6
    for kk=6.5:6.5
        for jj = 2.5:2.5
            plot([jj,jj-.5],[kk,ii],'k-'), hold on
        end
    end
end

for ii = 10:14
    for kk=6.5:6.5
        for jj = 2.5:2.5
            plot([jj,jj-.5],[kk,ii],'k-'), hold on
        end
    end
end

for ii = 1:6
    for kk=7.5:7.5
        for jj = 2.5:2.5
            plot([jj,jj-.5],[kk,ii],'k-'), hold on
        end
    end
end

for ii = 10:14
    for kk=7.5:7.5
        for jj = 2.5:2.5
            plot([jj,jj-.5],[kk,ii],'k-'), hold on
        end
    end
end

for ii = 1:6
    for kk=8.5:8.5
        for jj = 2.5:2.5
            plot([jj,jj-.5],[kk,ii],'k-'), hold on
        end
    end
end

for ii = 10:14
    for kk=8.5:8.5
        for jj = 2.5:2.5
            plot([jj,jj-.5],[kk,ii],'k-'), hold on
        end
    end
end

for ii = 1:6
    for kk=9.5:9.5
        for jj = 2.5:2.5
            plot([jj,jj-.5],[kk,ii],'k-'), hold on
        end
    end
end

for ii = 10:14
    for kk=9.5:9.5
        for jj = 2.5:2.5
            plot([jj,jj-.5],[kk,ii],'k-'), hold on
        end
    end
end

for ii = 1:6
    for kk=10.5:10.5
        for jj = 2.5:2.5
            plot([jj,jj-.5],[kk,ii],'k-'), hold on
        end
    end
end

for ii = 10:14
    for kk=10.5:10.5
        for jj = 2.5:2.5
            plot([jj,jj-.5],[kk,ii],'k-'), hold on
        end
    end
end

plot(2.5,6.5,'ko','markersize',8,'MarkerFaceColor','k')
plot(2.5,7.5,'ko','markersize',8,'MarkerFaceColor','k')
plot(2.5,8.5,'ko','markersize',8,'MarkerFaceColor','k')
plot(2.5,9.5,'ko','markersize',8,'MarkerFaceColor','k')
plot(2.5,10.5,'ko','markersize',8,'MarkerFaceColor','k')
plot(.5,6.5,'ko','markersize',8,'MarkerFaceColor','k')
plot(.5,7.5,'k.')
plot(.5,8.5,'k.')
plot(.5,9.5,'k.')
plot(.5,10.5,'ko','markersize',8,'MarkerFaceColor','k')

% plot([1,3],[])

ylim([-1 16])
xlim([0,3])
% axis off



