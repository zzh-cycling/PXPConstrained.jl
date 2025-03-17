using Yao

function toric_code_strings(m::Int, n::Int)
	li = LinearIndices((m, n))
	bottom(i, j) = li[mod1(i, m), mod1(j, n)] + m * n
	right(i, j) = li[mod1(i, m), mod1(j, n)]
	xstrings = Vector{Int}[]
	zstrings = Vector{Int}[]
	for i=1:m, j=1:n
		# face center
		push!(xstrings, [bottom(i, j-1), right(i, j), bottom(i, j), right(i-1, j)])
		# cross
		push!(zstrings, [right(i, j), bottom(i, j), right(i, j+1), bottom(i+1, j)])
	end
	return xstrings, zstrings
end

function toric_code_hamiltonian(m::Int, n::Int)
	xstrings, zstrings = toric_code_strings(m, n)
	sum([kron(2m*n, [i=>X for i in xs]...) for xs in xstrings[1:end-1]]) + sum([kron(2m*n, [i=>Z for i in zs]...) for zs in zstrings[1:end-1]])
end

h = toric_code_hamiltonian(3, 3)
