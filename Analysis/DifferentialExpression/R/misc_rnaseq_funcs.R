
plot_counts_matrix_log2_dist = function(matrix_file) {

	
	data = read.table(file=matrix_file, com='', row.names=1, header=T)

	conditions = colnames(data)
	colors = rainbow(length(conditions))


	plot(density(log2(data[,1])), col=colors[1], main=matrix_file, xlab='log2(frag_counts)', ylab='density')

	for (i in 2:length(data[1,])) {

		points(density(log2(data[,i])), type='l', col=colors[i])

	}

	legend('topright', conditions, col=colors, pch=15)

}

