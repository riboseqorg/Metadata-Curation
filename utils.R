termindexer <-  function(terms, dt) {
  out <- lapply(terms,  function(z) apply(dt,1, function(x) grep(z,x, ignore.case = T)))
  out1 <- lapply(out, function(x) which(0 != (x %>% lapply(length) %>% unlist)))
  out2 <- unlist(out1) %>% unique
  dt[out2]
}
