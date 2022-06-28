# Heli Juottonen 2019-2022 heli.juottonen@alumni.helsinki.fi
# R scripts for Elsi Hietaranta's willow manuscript
# https://github.com/helijuottonen/elsiwillow

# 1: PERMANOVA: testing for differences in bacterial and fungal community composition with treatments

# Disclaimer: I'm a self-taught R user. I try my very best to make sure that the scripts work, but I do not 
# guarantee that they are correct, make sense or follow any kind of best practices. Use at your own risk!

# RUN FIRST: 0_basic_processing.R

#library(vegan)

### permanovas for bagging, honeybees and site

# separating bag/no bag samples (because only no bag samples used for honeybee test)
#bacteria
elmetab.nb <- subset(elmetab, treat=="ctrl")
elotub.st.nb = elotub.st[row.names(elotub.st) %in% row.names(elmetab.nb),]
#fungi
elmetaf.nb <- subset(elmetaf, treat=="ctrl")
elotuf.st.nb = elotuf.st[row.names(elotuf.st) %in% row.names(elmetaf.nb),]

# defining site as a factor for later
elmetab$site <- as.factor(elmetab$site)
elmetaf$site <- as.factor(elmetaf$site)

# defining a variable for plant individuals within each site
elmetab$siterepl <- paste(elmetab$site, elmetab$repl, sep="_")
elmetaf$siterepl <- paste(elmetaf$site, elmetaf$repl, sep="_")

# setting permutations for treatment permanovas (within sites)
perm.b <- how(nperm = 999)
setBlocks(perm.b) <- with(elmetab, site)
perm.f <- how(nperm = 999)
setBlocks(perm.f) <- with(elmetaf, site)

# running permanovas (separate models for bagging (=treat), honeybees, site, plant individual)
# bacteria
adonis2(elotub.st ~ treat, distance="bray", permutations=perm.b, data=elmetab)
adonis2(elotub.st.nb ~ honeybees, distance="bray", permutations=999, data=elmetab.nb)
adonis2(elotub.st ~ site, distance="bray", permutations=999, data=elmetab)
adonis2(elotub.st ~ siterepl, distance="bray", permutations=perm.b, data=elmetab)

#adonis2(formula = elotub.st ~ siterepl, data = elmetab, permutations = perm.b, distance = "bray")
#Df SumOfSqs      R2      F Pr(>F)    
#siterepl 23   10.580 0.46614 2.6954  0.001 ***
#  Residual 71   12.117 0.53386                  
#Total    94   22.697 1.00000      

#fungi
adonis2(elotuf.st ~ treat, distance="bray", permutations=perm.f, data=elmetaf)
adonis2(elotuf.st.nb ~ honeybees, distance="bray", permutations=999, data=elmetaf.nb)
adonis2(elotuf.st ~ site, distance="bray", permutations=999, data=elmetaf)
adonis2(elotuf.st ~ siterepl, distance="bray", permutations=perm.f, data=elmetaf)

#adonis2(formula = elotuf.st ~ siterepl, data = elmetaf, permutations = perm.f, distance = "bray")
#Df SumOfSqs      R2     F Pr(>F)    
#siterepl 23   10.361 0.52177 3.368  0.001 ***
#  Residual 71    9.496 0.47823                 
#Total    94   19.857 1.00000   

### pollen permanova
adonis2(elotub.st.nb ~ pollen_remain, distance="bray", permutations=999, data=elmetab.nb)
adonis2(elotuf.st.nb ~ pollen_remain, distance="bray", permutations=999, data=elmetaf.nb)


### effect of plant individual (= 'repl' in metadata file)

# plant individuals are nested within sites 
# PERMANOVA in vegan can't handle nested designs -> splitting data by site

# splitting metadata of bacteria
elmetab.haa <- subset(elmetab, site=="Haasianlehto") # site A
elmetab.hp <- subset(elmetab, site=="Haapalehto") # site B
elmetab.juh <- subset(elmetab, site=="Juhanala") # cite C
elmetab.kek <- subset(elmetab, site=="Kekkila") # site D
elmetab.pal <- subset(elmetab, site=="Paloskyla") # site E
elmetab.pap <- subset(elmetab, site=="Pappinen") # site F

# splitting metadata of fungi
elmetaf.haa <- subset(elmetaf, site=="Haasianlehto") # site A
elmetaf.hp <- subset(elmetaf, site=="Haapalehto") # site B
elmetaf.juh <- subset(elmetaf, site=="Juhanala") # site C
elmetaf.kek <- subset(elmetaf, site=="Kekkila") # site D
elmetaf.pal <- subset(elmetaf, site=="Paloskyla") # site E
elmetaf.pap <- subset(elmetaf, site=="Pappinen") # site F

# splitting OTU data of bacteria
elotub.st.haa = elotub.st[row.names(elotub.st) %in% row.names(elmetab.haa),]
elotub.st.hp = elotub.st[row.names(elotub.st) %in% row.names(elmetab.hp),]
elotub.st.juh = elotub.st[row.names(elotub.st) %in% row.names(elmetab.juh),]
elotub.st.kek = elotub.st[row.names(elotub.st) %in% row.names(elmetab.kek),]
elotub.st.pal = elotub.st[row.names(elotub.st) %in% row.names(elmetab.pal),]
elotub.st.pap = elotub.st[row.names(elotub.st) %in% row.names(elmetab.pap),]

# splitting OTU data of fungi
elotuf.st.haa = elotuf.st[row.names(elotuf.st) %in% row.names(elmetaf.haa),]
elotuf.st.hp = elotuf.st[row.names(elotuf.st) %in% row.names(elmetaf.hp),]
elotuf.st.juh = elotuf.st[row.names(elotuf.st) %in% row.names(elmetaf.juh),]
elotuf.st.kek = elotuf.st[row.names(elotuf.st) %in% row.names(elmetaf.kek),]
elotuf.st.pal = elotuf.st[row.names(elotuf.st) %in% row.names(elmetaf.pal),]
elotuf.st.pap = elotuf.st[row.names(elotuf.st) %in% row.names(elmetaf.pap),]

# running permanovas on data split by site for bacteria
adonis2(elotub.st.haa ~ repl, distance="bray", permutations=999, data=elmetab.haa)
adonis2(elotub.st.hp ~ repl, distance="bray", permutations=999, data=elmetab.hp)
adonis2(elotub.st.juh ~ repl, distance="bray", permutations=999, data=elmetab.juh)
adonis2(elotub.st.kek ~ repl, distance="bray", permutations=999, data=elmetab.kek)
adonis2(elotub.st.pal ~ repl, distance="bray", permutations=999, data=elmetab.pal)
adonis2(elotub.st.pap ~ repl, distance="bray", permutations=999, data=elmetab.pap)

# running permanovas on data split by site for fungi
adonis2(elotuf.st.haa ~ repl, distance="bray", permutations=999, data=elmetaf.haa)
adonis2(elotuf.st.hp ~ repl, distance="bray", permutations=999, data=elmetaf.hp)
adonis2(elotuf.st.juh ~ repl, distance="bray", permutations=999, data=elmetaf.juh)
adonis2(elotuf.st.kek ~ repl, distance="bray", permutations=999, data=elmetaf.kek)
adonis2(elotuf.st.pal ~ repl, distance="bray", permutations=999, data=elmetaf.pal)
adonis2(elotuf.st.pap ~ repl, distance="bray", permutations=999, data=elmetaf.pap)

