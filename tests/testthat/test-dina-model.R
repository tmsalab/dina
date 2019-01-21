context("test-dina-model")

test_that("Verify Bad Input detected", {
    
    expect_error(dina(1:5, sim_q_matrix(10, 3)), info = "Detected Bad Item Data Matrix")
    expect_error(dina(sim_q_matrix(10, 3), 1:5), info = "Detected Bad Q Matrix")

})

test_that("Verify Model Reproducibility", {
    
    # Setup Parameters
    N = 15   # Number of Examinees / Subjects
    J = 10   # Number of Items
    K = 2    # Number of Skills / Attributes
    
    # Assign slipping and guessing values for each item
    ss = gs = rep(.2, J)
    
    # Simulate identifiable Q matrix
    Q = sim_q_matrix(J, K)
    
    # Simulate subject attributes
    subject_alphas = sim_subject_attributes(N, K)
    
    # Item data
    items_dina = sim_dina_items(subject_alphas, Q, ss, gs)
    
    ## Estimate two models ----
    # Set a seed for reproducibility
    set.seed(192)
    model1 = dina(items_dina, Q, chain_length = 100)
    
    # Set the same seed
    set.seed(192)
    model2 = dina(items_dina, Q, chain_length = 100)
    
    expect_equal(model1, model2,
                 info = "Models are reproducible under same seed.")
    
})

