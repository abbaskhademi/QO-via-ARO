%  Liuzzi et al. (2019), A new branch-and-bound algorithm for standard quadratic programming problems

function myrandom(r,a,b)
	r = mod(r*41475557.,1.)
	return (r,a+r*(b-a))
end

n      = 99
dens   = 0.9
dvert  = 2
myseed = 1

for igen = 0:5
	r = (4.*myseed + 1.)/16384./16384.
	myseed = myseed+1

	Fc = zeros(n+1,n+1)
	Fl = zeros(n+1,n+1)
	F  = zeros(n+1,n+1)

	for i = 1:n
		for j = (i+1):(n+1)
			(r,num) = myrandom(r,0,1)
			if num < dens
				(r,num) = myrandom(r,0.,10.)
				Fc[i,j] = Fc[j,i] = num
			else
				(r,num) = myrandom(r,-10.,0.)
				Fc[i,j] = Fc[j,i] = num
			end
		end
	end
	for i = 1:(n+1)
		(r,num) = myrandom(r,0,dvert)
		Fl[i,i] = num
	end
	for i = 1:n
		for j = (i+1):(n+1)
			Fl[i,j] = Fl[j,i] = 0.5*(Fl[i,i] + Fl[j,j])
		end
	end
	for i = 1:(n+1)
		for j = i:(n+1)
			F[i,j] = F[j,i] = Fl[i,j] - Fc[i,j]
		end
	end
	
	writedlm("Problem_"*string(n+1)*"x"*string(n+1)*"("*string(dens)*")_"*string(igen)*".txt",F)
end

nothing
