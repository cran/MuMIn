
.progbar <-
function(min = 0, max = 1, initial = 0, title = "")
    local({
       if(.Platform$GUI == "Rgui") {
			pb <- utils::winProgressBar(min = min, max = max,, title = title)
            setpb <- utils::setWinProgressBar
        } else {
            pb <- utils::txtProgressBar(min = min, max = max, style = 3L)
            setpb <- utils::setTxtProgressBar
        }
        function(...) setpb(pb, ...)
        #function(value) koCmd(sprintf("kor.setProgressBar(%d, %f, %s, %s)", 
         #   pb$pb, as.double(value), "''", "''"))
    })
    
.closeprogbar <- function(x) close(environment(x)$pb)
