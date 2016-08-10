library(edgeR)



### 3.3.2 Nested interaction formulas

targets <- data.frame(Treat = rep(c("Placebo", "Drug"), each = 6), Time = c(rep(c("0h", "1h", "2h"), each = 2), rep(c("0h", "1h", "2h"), each = 2)))

rownames(targets) <- paste(1:nrow(targets), targets$Treat, targets$Time, sep = ".")


targets$Treat <- relevel(targets$Treat, ref="Placebo")


model.matrix(~ Treat + Treat:Time, data=targets)

model.matrix(~ Treat + Time + Treat:Time, data=targets)




model.matrix(~ Treat + Time + Time:Treat, data=targets)

model.matrix(~ Time + Treat + Time:Treat, data=targets)

model.matrix(~ Treat + Time:Treat, data=targets)













