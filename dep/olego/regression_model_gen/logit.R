args = commandArgs();
input_file=args[4];
output_file=args[5];
sink(output_file);
data=read.csv(input_file);
summary(data);
attach(data);
#names(data);
fit = glm(class~size+score, family=binomial(link="logit"))
summary(fit)


#if ( file.info(input_file)["size"]>0 )
#{
#data=read.csv("logit.csv");
#attach(data);
#fit = glm(class~size+score, family=binomial(link="logit"))
#p[i]=dbeta(i/1000,f[["estimate"]][["shape1"]],f[["estimate"]][["shape2"]])

#write(p,file=paste(input_file,"fit",sep="."),ncolumns=1);

