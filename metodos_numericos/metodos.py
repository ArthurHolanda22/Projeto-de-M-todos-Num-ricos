def bissecao(func, a, b, tol=1e-6, max_iter=100):
    if func(a) * func(b) >= 0:
        raise ValueError("O intervalo [a, b] não contém troca de sinal válida.")
        
    iteracoes = 0
    historico_erro = []
    
    for _ in range(max_iter):
        iteracoes += 1
        c = (a + b) / 2.0
        fc = func(c)
        erro = (b - a) / 2.0
        historico_erro.append(erro)
        
        if abs(fc) < tol or erro < tol:
            return c, iteracoes, historico_erro
            
        if func(c) * func(a) < 0:
            b = c
        else:
            a = c
            
    return (a + b) / 2.0, iteracoes, historico_erro


def newton_raphson(func, dfunc, x0, tol=1e-6, max_iter=100):
    x = x0
    historico_erro = []
    
    for iteracoes in range(1, max_iter + 1):
        fx = func(x)
        dfx = dfunc(x)
        
        if abs(dfx) < 1e-12:
            raise ZeroDivisionError("Derivada próxima de zero.")
            
        x_novo = x - fx / dfx
        erro = abs(x_novo - x)
        historico_erro.append(erro)
        
        if erro < tol or abs(fx) < tol:
            return x_novo, iteracoes, historico_erro
            
        x = x_novo
        
    return x, max_iter, historico_erro


def secante(func, x0, x1, tol=1e-6, max_iter=100):
    historico_erro = []
    
    for iteracoes in range(1, max_iter + 1):
        f0, f1 = func(x0), func(x1)
        
        if abs(f1 - f0) < 1e-12:
            raise ZeroDivisionError("Diferença f(x1) - f(x0) muito pequena.")
            
        x_novo = x1 - f1 * (x1 - x0) / (f1 - f0)
        erro = abs(x_novo - x1)
        historico_erro.append(erro)
        
        if erro < tol or abs(f1) < tol:
            return x_novo, iteracoes, historico_erro
            
        x0, x1 = x1, x_novo
        
    return x1, max_iter, historico_erro