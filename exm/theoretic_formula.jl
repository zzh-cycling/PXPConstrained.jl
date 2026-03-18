EBA(t)=[cos(t/2)^2 0 0 cos(t/2)^2*sin(t/2)^2; 0 0 0 0; 0 0 0 0; sin(t/2)^4 0 0 sin(t/2)^4]
EAB(t)=[sin(t/2)^2 0 0 cos(t/2)^2*sin(t/2)^2; 0 0 0 0; 0 0 0 0; cos(t/2)^4 0 0 cos(t/2)^4]

EOBAelement(t) = -2cos(t/2)^3*sin(t/2)
EOABelement(t) = 2sin(t/2)^3*cos(t/2)
EOBA(t) = [EOBAelement(t) 0 0 0;0 0 0 0;0 0 0 0; 0 0 0 0]
EOAB(t) = [EOABelement(t) 0 0 0;0 0 0 0;0 0 0 0; 0 0 0 0]

f(t) = sqrt(2)*sqrt(44*cos(2*t)-3*cos(4*t)+87)
λ1(t) = 1/32*(14+2cos(2*t)-f(t))
λ2(t) = 1/32*(14+2cos(2*t)+f(t))

vec1RBA(t) = [1/32*(-32+14csc(t/2)^4+2cos(2*t)*csc(t/2)^4-csc(t/2)^4*f(t)),0,0,1]
vec2RBA(t) = [1/32*(-32+14csc(t/2)^4+2cos(2*t)*csc(t/2)^4+csc(t/2)^4*f(t)),0,0,1]
vec1LBA(t) = [1/32 * (-32 + 14 * csc(t/2)^4 + 2 * cos(2 * t) * csc(t/2)^4 - f(t) * csc(t/2)^4 )*tan(t/2)^2,0,0,1]
vec2LBA(t) = [1/32 * (- 32 + 14 * csc(t/2)^4 + 2 * cos(2 * t) * csc(t/2)^4 + f(t) * csc(t/2)^4 )*tan(t/2)^2 ,0,0,1]

vec1RAB(t) = [-1/32*(-14+32cos(t/2)^4-2cos(2*t)+f(t))*sec(t/2)^4,0,0,1]
vec2RAB(t) = [1/32*(-32+14sec(t/2)^4+2cos(2*t)*sec(t/2)^4+sec(t/2)^4*f(t)),0,0,1]
vec1LAB(t) = [-1/32*(-14+32cos(t/2)^4-2cos(2*t)+f(t))*csc(t/2)^2*sec(t/2)^2,0,0,1]
vec2LAB(t) = [-1/32*(-14+32cos(t/2)^4-2cos(2*t)-f(t))*csc(t/2)^2*sec(t/2)^2,0,0,1]

nmvec(st::Vector{Float64}) = st ./ norm(st)


normvec1RBA(t) = nmvec(vec1RBA(t))
normvec2RBA(t) = nmvec(vec2RBA(t))
normvec1LBA(t) = vec1LBA(t)./(vec1LBA(t)'*normvec1RBA(t))
normvec2LBA(t) = vec2LBA(t)./(vec2LBA(t)'*normvec2RBA(t))

normvec1RAB(t) = nmvec(vec1RAB(t))
normvec2RAB(t) = nmvec(vec2RAB(t))
normvec1LAB(t) = vec1LAB(t)./(vec1LAB(t)'*normvec1RAB(t))
normvec2LAB(t) = vec2LAB(t)./(vec2LAB(t)'*normvec2RAB(t))

ρ3(t) = 32*4/((14+2*cos(2*t)+f(t))*f(t)^2*csc(t)^6)