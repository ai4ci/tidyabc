browse_on_error = function(expr, env = rlang::caller_env()) {
  expr = rlang::enexpr(expr)
  tmp = try(eval(expr, envir = env), silent = TRUE)
  if (inherits(tmp, "try-error")) {
    message("Code block threw: ", attr(tmp, "condition")$message)
    message("Entering debug: ")
    expr = c(expression(browser()), expr)
    eval(expr, envir = env)
  }
  return(tmp)
}
