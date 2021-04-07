#=
	test hclust to newick tree format

=#

using Clustering
include("../src/hclust2newick.jl")

@testset "hclust2newick()" begin
    A = [ 
		0 0; 1 2; 1 1; 1 0; 2 1; 9 1; 10 2; 10 1; 10 0; 11 1; 13 1; 
		14 2; 14 1; 14 0; 15 1; 11 4; 12 5; 12 4; 12 3; 13 4
	]
	labels = [
		"B", "B", "B", "B", "B", "M", "M", "M", "M", "M", 
		"P", "P", "P", "P", "P", "G", "G", "G", "G", "G"
	]
	n = first(size(A))
	D = zeros(n, n)

	for i in 1:n-1
		for j in i+1:n
			D[i,j] = D[j,i] = D[i,j] = sqrt(sum((A[i,:] .- A[j,:]).^2))
		end 
	end
    @test hclust2newick(hclust(D); labels=labels) == "((B:0.0,(B:0.0,(B:0.0,(B:0.5,B:0.5):0.0):0.0):0.0):0.0,((G:0.0,(G:0.0,(G:0.0,(G:0.5,G:0.5):0.0):0.0):0.0):0.0,((P:0.0,(P:0.0,(P:0.0,(P:0.5,P:0.5):0.0):0.0):0.0):0.0,(M:0.0,(M:0.0,(M:0.0,(M:0.5,M:0.5):0.0):0.0):0.0):0.0):0.0):0.0);"

end



