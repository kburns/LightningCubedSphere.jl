"reciprocal_log(z, c, s, a) = ..."
function reciprocal_log(z, c, s, a)
    f = zero(z)
    for (aj, sj) in zip(a, s)
        for k = 0:3
            # this step symmetrizes the branch cut
             fjk = @. im^k * cispi(1/4) * aj / (log(-(z/(im^k * c) - 1)) + im*pi - sj)
             f += fjk
             fjk = @. im^k * cispi(1/4) * aj / (log(-(z/(im^k * c) - 1)) - im*pi - sj)
             f += fjk
         end
    end
    return f
end

"Log-lightning representation."
log_lightning(z, c, s, a, b) = reciprocal_log(z, c, s, a) + runge(z, b)

# Log least squares fit
if include_log
    s = LinRange(-5,20,na)
    Re_fLog(q) = real(log_lightning(z, c, s, q[1:na], q[na+1:end]))
    qLog = fit(Re_fLog, ones(ns), zeros(na+nb))
    println(" (fLog)")
    aLog = qLog[1:na]
    bLog = qLog[na+1:end]
end

if include_log
    fLog_r = log_lightning(z_r, c, s, aLog, bLog)
    @printf "resampling error: %.2e (fLog)\n" maximum(abs.(real(fLog_r) .- 1))
end
