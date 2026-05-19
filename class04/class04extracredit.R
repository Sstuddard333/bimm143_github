source("http://thegrantlab.org/misc/cdc.R")
View(cdc)
head(cdc$height, 7)
tail(cdc$weight, 20)
plot(cdc$height, cdc$weight, xlab = "Height", ylab = "Weight")
cor(cdc$height, cdc$weight)
weight_kg <- cdc$weight * 0.454
height_m <- cdc$height * 0.0254

BMI <- (weight_kg) / (height_m^2)
plot (cdc$height, BMI, xlab = "Height", ylab = "BMI")
cor(cdc$height, BMI)


sum(BMI >= 30)
plot(cdc[1:100,5], cdc[1:100,6])

obeseBMI <- BMI >= 30
table( cdc$gender, obeseBMI)
