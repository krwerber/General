attach(housingdata)

skim(sale_price)
summary(sale_price)

housingdata %<>%
  filter(!is.na(sale_price))


housingdata$sale_price %<>%
  parse_number()

n <- nrow(housingdata)

skim(housingdata)

housingdata[, which(colSums(is.na.data.frame(housingdata)) == n)] <- NULL

housingdata$is_coop <- coop_condo == "co-op"

cats_allowed <- cats_allowed == "yes"

dogs_allowed <- dogs_allowed == "yes"

housingdata$pets_allowed <- cats_allowed | dogs_allowed

housingdata$cats_allowed <- NULL
housingdata$dogs_allowed <- NULL

colnames(housingdata)

skim(housingdata)

housingdata[, 1:12] <-  NULL
housingdata[, 2:11] <- NULL

housingdata$zip_code <- stringr::str_extract(housingdata$full_address_or_zip_code,
                                                        pattern = "11[0-9]{3}") %>%
                                            as.factor()

housingdata$full_address_or_zip_code <- NULL

garage_exists <-  !is.na(garage_exists)
housingdata$garage_exists <- garage_exists

housingdata$date_of_sale %<>%
  as.Date(format = "%m/%d/%Y")

housingdata$coop_condo <- NULL

housingdata$model_type <- NULL

housingdata$num_floors_in_building <- NULL

housingdata$pct_tax_deductibl <- NULL

housingdata$common_charges <- NULL

apply(housingdata$dining_room_type, FUN = function(x){if(x == "dining area"){x <- "other"}})

housingdata$dining_room_type[housingdata$dining_room_type == "dining area"] <- "other"

housingdata$fuel_type[housingdata$fuel_type %in% c("none", "other", "Other")] <- "other"

housingdata$kitchen_type[kitchen_type %in% c("eat in", "Eat in", "Eat In")] <- "eat in"

housingdata$kitchen_type[kitchen_type %in% c("combo", "Combo")] <- "combo"

housingdata$kitchen_type[kitchen_type == "1955"] <- NA

housingdata$kitchen_type <- droplevels.factor(x = housingdata$kitchen_type)

housingdata$num_half_bathrooms[is.na(num_half_bathrooms)] <- 0
housingdata$parking_charges[is.na(parking_charges)] <- 0

housingdata$maintenance_cost %<>%
  parse_number()

housingdata$parking_charges %<>%
  parse_number()

housingdata$common_charges %<>%
  parse_number()

housingdata$total_taxes %<>%
  parse_number()

housingdata$dining_room_type %<>%
  as.factor()

housingdata$WorkerId %<>%
  as.factor()

housingdata$community_district_num %<>%
  as.factor()

housingdata$fuel_type %<>%
  as.factor()

housingdata$kitchen_type %<>%
  as.factor()

housingdata$zip_code %<>%
  as.factor()

housingdata$community_district_num %<>%
  as.factor()

housingdata$approx_year_built %<>%
  as.integer()

ggplot(data = housingdata, aes(approx_year_built, sale_price)) +
  geom_point()

ggplot(data = housingdata, aes(date_of_sale, sale_price)) +
  geom_point()

housingdata$date_of_sale <- NULL

M <- tbl_df(apply(is.na(housingdata), 2, as.numeric))
colnames(M) = paste("is_missing_", colnames(housingdata), sep = "")

M <- tbl_df(t(unique(t(M))))
M %<>% 
  select_if(function(x){sum(x) > 0})

Xnum <- housingdata %>%
          select_if(function(x){is.logical(x) | is.numeric(x)}) %>%
          select(-sale_price) %>% 
          as.matrix()
  
Xnum_imp <- missForest::missForest(Xnum)
  
housing_imp_1 <- cbind.data.frame(Xnum_imp$ximp, housingdata %>%
                     select_if(function(x){!(is.logical(x) | is.numeric(x))}))


Xcat <- housingdata %>%
          select_if(is.factor)

housing_imp_final <- mice::mice(housing_imp_1, method = "polyreg")

X <- complete(housing_imp_final, sample(1:5, size = 1)) %>%
  cbind(M)

D <- cbind(X, sale_price)

X_reg <- D %>%
  select(-WorkerId, -community_district_num, -zip_code, -sale_price) 

n <- nrow(D)

train_indices <- sample(1:n, size = 0.8 * n, replace = FALSE)
test_indices <- setdiff(1:n, train_indices)

testthat::expect_equal(sort(c(train_indices, test_indices)), 1:n)

X_train <- X_reg[train_indices, ]
y_train <- sale_price[train_indices]
X_test <- X_reg[test_indices, ]
y_test <- sale_price[test_indices]

D_train <- cbind(X_train, y_train)
D_test <- cbind(X_test, y_test)

reg <- lm(y_train ~ ., data = D_train)
summary(reg)

lin_oos_rmse <- sqrt(mean((y_test - predict(reg, X_test))^2))
oos_rsq <- 1 - sum((y_test - predict(reg, X_test))^2) / sum((y_test - mean(y_test))^2)

options(java.parameters = "-Xmx4000m")
tree_reg <- YARFCART(X_train, 
                     y_train, 
                     bootstrap_indices = 1:nrow(X_train))

part_reg <- rpart(sale_price ~ ., data = D, method = "anova")

printcp(part_reg)
plot(part_reg, uniform = TRUE)
text(part_reg, use.n=TRUE, all=TRUE, cex=.8)

illustrate_trees(tree_reg, max_depth = 4)

tree_rmse <- sqrt(mean((y_test - predict(tree_reg, X_test))^2))

rf_mod <- YARF(X_train, y_train, num_trees = 1000)

print(c("Random Forest RMSE:", rf_mod$rmse_oob))

tree_reg$rmse_oob
rf_resids <- cbind.data.frame(y_test, resids <- y_test - predict(rf_mod, X_test))

ggplot(data = rf_resids, aes(y_test, resids)) +
  geom_point()


