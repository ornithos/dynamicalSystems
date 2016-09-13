assert <- function(condition, reason, ...) {
    
    if(!condition) {
        stop(sprintf(paste0("\n",reason), ...))
    }
}
        