@testset "UProblem" begin
    
    #Example QP
    x_0 = ones(3);
    f(x) = P.dot(x, x);
    function g!(gx, x)
        @. gx = 2 * x;
        return nothing
    end
    up = UProblem(x_0, f, g!);
    @test up.x_0 == x_0;
    @test up.f(up.x_0) == 3.0;
    gx_0 = zeros(3);
    up.g!(gx_0, up.x_0);
    @test gx_0 == 2.0 * ones(3);

end