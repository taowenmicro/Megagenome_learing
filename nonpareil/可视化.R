
source("./Nonpareil.R")
png("./nonpareil_forward.png")

# install.packages("Nonpareil")
library(Nonpareil)
np <- Nonpareil.curve("../nonpareil/nonpareil_forward.npo")
legend('bottomright', legend = c(paste("Coverage: ", round(np$C*100,digits=2), "%"),paste("Actual effort =", round(np$LR/1000000,digits=2), "Mbp"), paste("Required effort for 95% coverage=", round(np$LRstar/1000000,digits=2), "Mbp"), paste("Diversity =", round(np$diversity,digits=2))),bty = "n") 
dev.off()
