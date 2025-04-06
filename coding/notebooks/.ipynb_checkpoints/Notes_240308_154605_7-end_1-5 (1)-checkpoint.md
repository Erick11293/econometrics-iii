(2) To get $y_{i}$ as a function of the past values, we need to invert the MA process If we define $\gamma=-\gamma_{1}$

$$
\begin{gathered}
(1-\phi L) y_{i}=(1-\gamma L) \varepsilon_{t} \\
\frac{(1-\phi L)}{(1-\gamma l)} y_{t}=\varepsilon_{t} \\
(1-\phi L)_{i=0}^{\infty}(\gamma L)^{i} y_{t}=\varepsilon_{t} \\
\left.(1-\phi L) \sum_{i=0}^{t-1}(\partial L)^{i} y_{t}+\sum_{i=t}^{\infty}(\gamma l)^{i} y_{t}\right]=\varepsilon_{t} \\
(1-\phi L) \sum_{i=0}^{t-1}(\gamma l)^{i} y_{t}+\gamma^{t} \varepsilon_{0}=\varepsilon_{\tau} \\
y_{t}+\sum_{i=1}^{t-1}(\gamma L \geqslant 2)^{i} y_{t}-\sum_{i=0}^{t-1} \phi \gamma^{i}(l)^{i+1} y_{t}+\gamma^{t} \varepsilon_{0}=\varepsilon_{t}
\end{gathered}
$$

$$
\begin{aligned}
& y_{t}+\sum_{i=1}^{t-1}(\partial L)^{i} y_{t}-\sum_{t=1}^{t-1} \phi \gamma^{i-1} L^{i} y_{t} \\
&+\gamma^{t} \varepsilon_{0}-\phi \gamma^{t-1} y_{0}=\varepsilon_{t} \\
& y_{t}-\sum_{i=1}^{t-1} \gamma^{i-1}(\phi-\gamma) L^{i} y_{t}=\varepsilon_{t}-\gamma^{t} \varepsilon_{0} \\
&+\phi \gamma^{t-1} y_{0} \\
& y_{t}=\sum_{i=1}^{t-1} \gamma^{i-1}(\phi-\gamma) L^{i} y_{t}+\varepsilon_{t} \\
&-\gamma^{t} \varepsilon_{0}+\phi \gamma^{t-1} y_{0} \\
& x_{i} \geqslant 2 \\
& y_{1}=\varepsilon_{1}-\gamma \varepsilon_{0}+\phi y_{0} \quad i=1
\end{aligned}
$$

Then, the expectation for $k \geqslant 1$ is:

$$
E\left(y_{t+k} I^{t}, \varepsilon_{0}\right)=\phi^{k-1} E\left(y_{t+1} \mid I^{t}, \varepsilon_{0}\right)
$$

$\forall K \geqslant 1$, as $\quad\left(\varepsilon_{t+K}\left(I^{\prime}, \varepsilon_{0}\right)=0\right.$. So:

$$
E\left(y_{t+k} \mid I_{1}^{t}, \varepsilon_{0}\right)=\phi^{k-1}\left(\sum_{i=0}^{t-1} \partial^{i}(\phi-\partial) y_{t-1}+\partial^{t+1} \varepsilon_{0}+\phi \gamma^{t} y_{0}\right)
$$

(3) The conditional login Kelli hood would be normal. We know the conch tuinal expectation of $y_{t+1} \mid I^{t}$, so we can use it to write the likelihood

$$
y_{t+1} \mid I^{t} \sim N\left(E\left(y_{t+1} \mid I^{t}, \varepsilon_{0}\right), \sigma^{2}\right)
$$

So

$$
P\left(\left\{y_{t}\right\}_{t=1}^{T}\right)=\prod_{i=1}^{T} \hat{P}\left(y_{t} \mid I^{t-1}, \varepsilon_{0}\right)
$$

