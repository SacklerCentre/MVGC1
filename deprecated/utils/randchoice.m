function co = randchoice(things,c)

perm = randperm(length(things));
co = things(perm(1:c));
