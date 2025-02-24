#---------------------#
#------- TESTS -------#
#---------------------#
test_that("update_velocity does not change when V is zero", {

  N <- 5
  d <- 3
  X <- matrix(1, nrow = N, ncol = d)  # all positions are the same
  # TEST setup: all particles are not moving
  V <- matrix(0, nrow = N, ncol = d)
  P <- X  # personal bests are the same as current positions
  informants <- lapply(1:N, function(i) 1:N)

  # perform the velocity update
  new_V <- update_velocity(X, V, P, informants, w, c)

  # test if the velocities remain unchanged
  expect_equal(new_V, V, info = "update_velocity does not change when V is zero.")
})
test_that("non-randomness in velocity change proportional to `w `when position=best", {

  # test setup: all particles have identical positions and fitness values
  N <- 5
  d <- 3
  X <- matrix(1, nrow = N, ncol = d)  # all positions are the same
  V <- matrix(1, nrow = N, ncol = d)  # initial velocities
  P <- X  # personal bests are the same as current positions
  informants <- lapply(1:N, function(i) 1:N)
  w = 3
  # perform the velocity update
  new_V <- update_velocity(X, V, P, informants, w, c)
  # by the algebra... \forall j \in \mathcal{I}_i, (\mathbf{p}_t^{(j)} - \mathbf{X}_t^{(i)}) = 0
  # hence,
  # V_{t + 1} = w * V_t
  expect_equal(new_V, w*V, info = "on-randomness in velocity change proportional to `w `when position=best.")
})
