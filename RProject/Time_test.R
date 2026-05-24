
# Number of repetitions
repetitions <- 1

# Vector to store elapsed times
elapsed_times <- numeric(repetitions)

# Run the function multiple times
for (i in 1:repetitions) {
  start_time <- Sys.time()
  Omega_statistic(woodmouse[, 1:643])
  end_time <- Sys.time()
  elapsed_times[i] <- end_time - start_time
}

# Calculate the average elapsed time
average_elapsed_time <- mean(elapsed_times)

# Print the average elapsed time
print(average_elapsed_time)

rm(elapsed_times, average_elapsed_time)




# Start timing
start_time <- Sys.time()

# Run the function
Omega_statistic(woodmouse[, 1:643])

# End timing
end_time <- Sys.time()

# Calculate the elapsed time
elapsed_time <- end_time - start_time

# Print the elapsed time
print(elapsed_time)





# Create the data frame with 4 columns and 3 rows
df <- data.frame(
  Column1 = c("Function", "Achaz_correction", "Omega_statistic"),
  Column2 = c(322, 0.00346365, 34.17451),
  Column3 = c(643, 0.00428102, 69.906),
  Column4 = c(965, 0.009723878, 107.097)
)

# Set the first row as column names
names(df) <- df[1, ]

# Remove the first row
df <- df[-1, ]

# Set the first column as row names
row.names(df) <- df[, 1]

# Remove the first column
df <- df[, -1]

# Print the data frame
print(df)


ggplot(df, aes(x = variable, y = value, color = RowNames)) +
  geom_line() 





# Convert tibble to data frame
Book <- as.data.frame(Book)

# Set the first column as row names
rownames(Book) <- Book[, 1]

# Remove the first column
Book <- Book[, -1]




Book_long <- tidyr::gather(Book, key = "Category", value = "Value")

# Plot using ggplot2
ggplot(Book_long, aes(x = names(Book_long)[3], y = Value, color = Category)) +
  geom_line() +
  labs(x = "X Axis Label", y = "Y Axis Label", title = "Title") +
  theme_minimal()
