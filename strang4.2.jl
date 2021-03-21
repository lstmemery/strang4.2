### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 113d38f4-899a-11eb-27ef-bd783b2a1ba2
using LinearAlgebra

# ╔═╡ dbeaa654-8a55-11eb-21f7-477091d9849b
x_hat(A, b) = (A'*A)\A'*b

# ╔═╡ e71e2b92-8a56-11eb-3fe1-09ec3fe2771c
p(A, b) = A * x_hat(A, b)

# ╔═╡ f5f02be8-8a56-11eb-1228-35ab779e3be2
e(A, b) = b - p(A, b)

# ╔═╡ 0fb707e0-8a57-11eb-250a-959e36596fdb
P(A) = A* inv((A' * A)) * A'

# ╔═╡ 6a461b84-8a56-11eb-2a35-912dfb4ca40d
A = [
	1 0;
	1 1;
	1 2
]

# ╔═╡ 7ec7c1fc-8a56-11eb-2866-c39711f9787e
b = [6 0 0]'

# ╔═╡ 832c7ac6-8a56-11eb-245e-fdb85294c95b
x_hat(A, b)

# ╔═╡ b6f7a2f4-8a56-11eb-0d45-29ecacfd5b93
p(A,b)

# ╔═╡ 2296918c-8a57-11eb-25ca-7f1d9b3a6d6f
e(A,b)

# ╔═╡ 93e57e9e-8a56-11eb-0dd8-cd087d26a743
P(A) *6

# ╔═╡ 2bd4456e-8a57-11eb-0d27-87a4529d9b3c
md"""
## Q1
"""

# ╔═╡ 5f38e6bc-8a57-11eb-11ff-d101b514d5f5
q1aa = [1 1 1]'

# ╔═╡ 6e407d1e-8a57-11eb-3da2-0b5ca81c87a1
q1ab = [1 2 2]'

# ╔═╡ 754de3f8-8a57-11eb-0871-31b59419fe80
p(q1aa, q1ab)

# ╔═╡ 9030867c-8a59-11eb-0dd4-233a23acba24


# ╔═╡ 8541822c-8a59-11eb-17eb-5936e9f7388d


# ╔═╡ 7d481752-8a57-11eb-33ed-dd651cc0f199
q1ae = e(q1aa, q1ab)

# ╔═╡ 8a2fa55c-8a57-11eb-2ac7-3d295325c893
q1ae' * q1aa

# ╔═╡ a3971ff4-8a57-11eb-2382-a1da4802eb2c
q1ba = [-1 -3 -1]'

# ╔═╡ 17325ac8-8a58-11eb-3963-69200384e81f
q1bb = [1 3 1]'

# ╔═╡ 1ee9963a-8a58-11eb-338a-cf27c7fe7d74
p(q1ba, q1bb)

# ╔═╡ 2fabe29a-8a58-11eb-0cb9-950bf424bda0
q1be = e(q1ba, q1bb)
# any vector multiplied by the 0 vector will be zero

# ╔═╡ 4335141e-8a58-11eb-3f07-f74ccfd631d9
md"""
## Q2
"""

# ╔═╡ 9144618a-8a58-11eb-0664-bb1f527b3a93
q2aa = [1 0]'

# ╔═╡ f5696b14-8a59-11eb-1185-57cecbd9e437
q2aP = P(q2aa)

# ╔═╡ 33989bbc-8a5a-11eb-0d8e-39c297ea0627
q2ab = [cos(pi) sin(pi)]' # Selecting pi/2 arbitrarily

# ╔═╡ 877bb944-8a5a-11eb-1cf8-414885322af4
q2aP * q2ab

# ╔═╡ 9281c784-8a5a-11eb-3cb7-179e43897555
q2ba = [1 -1]'

# ╔═╡ b13d1144-8a5a-11eb-113b-c15cd1ec4a65
q2bb = [1 1]'

# ╔═╡ b97a79c6-8a5a-11eb-0009-a3b0a3dc6b45
p(q2ba, q2bb)

# ╔═╡ c3c2c33e-8a5a-11eb-2d7f-47784785ff89
md"""
## Q3
"""

# ╔═╡ f69a69b0-8a5a-11eb-25d9-6374163cd317
q3aP = P(q1aa)

# ╔═╡ 0d172a2a-8a5b-11eb-35bd-958247b6a94c
q3aP * q3aP

# ╔═╡ 18b99b60-8a5b-11eb-35c7-a30bf4517452
q3aP * q1ab

# ╔═╡ 2651f81c-8a5b-11eb-07b0-df1b53c0f16f
q3bP = P(q1ba)

# ╔═╡ 5d4e180a-8a5b-11eb-3658-bb32f621f04b
q3bP * q3bP

# ╔═╡ 67ee890c-8a5b-11eb-2845-51dfba63507c
q3bP * q1bb

# ╔═╡ 827e5ad6-8a5b-11eb-13da-cdfddadb87da
md"""
## Q4

(I already did a!)
"""

# ╔═╡ 803ddbbe-8a5b-11eb-18da-b7da029a7fb0
q4bP = P(q2ba)

# ╔═╡ 0eed69da-8a5c-11eb-3313-052ddd989a8f
q2aP + q4bP

# ╔═╡ f3edfa3c-8a5b-11eb-23a1-0ff1b2a8c811
(q2aP + q4bP)^2 # It is not true.

# ╔═╡ 073eb64e-8a5c-11eb-1609-4799913b89ff
md"""
## Q5
"""

# ╔═╡ 767e04de-8a5b-11eb-3645-2d2ce2984b8c
q5a1 = [-1 2 2]'

# ╔═╡ 44b4218a-8a5c-11eb-212f-9f128a8ec531
q5a1P = P(q5a1)

# ╔═╡ 5a9d18ce-8a5b-11eb-28c7-014c17246015
q5a2 = [2 2 -1]'

# ╔═╡ 6b5683b4-8a5c-11eb-1410-13b41f925244
q5a2P = P(q5a2)

# ╔═╡ 515a7384-8a5b-11eb-1ee9-8d10c0d2a66e
q5a1P*q5a2P
# This matrix is 0 because...

# ╔═╡ 8c2271a2-8a5c-11eb-3e47-53ce1261473e
q5a1' * q5a2 # The vectors are orthogonal!

# ╔═╡ a1f04194-8a5c-11eb-33fc-a3ce4730ce7f
md"""
## Q6
"""

# ╔═╡ b7d0ad6e-8a5c-11eb-0c90-0197079e8e9c
q6b = [1 0 0]'

# ╔═╡ c36aa0ee-8a5c-11eb-0806-d3f0366148eb
q6a = [2 -1 2]'

# ╔═╡ ce18fc02-8a5c-11eb-004d-e11309767ac4
q6P = P(q6a)

# ╔═╡ dd96491e-8a5c-11eb-24a6-2f4c230c0311
q6p1 = q5a1P * q6b

# ╔═╡ f45ac2aa-8a5d-11eb-2913-593cf0569cd3
q6p2 = q5a2P * q6b

# ╔═╡ 04281bc4-8a5e-11eb-2cd9-774a44c2a49f
q6p3 = q6P * q6b

# ╔═╡ ed83a4b0-8a5d-11eb-00c1-d575a68861ba
q6p1 + q6p2 + q6p3 # This is equal to b

# ╔═╡ eac9248e-8a5d-11eb-35cb-8585b1265fb5
md"""
## Q7
"""

# ╔═╡ 4c449b28-8a5e-11eb-36f0-77625dc02ded
q5a1P + q5a2P + q6P

# ╔═╡ 7149ec0c-8a5e-11eb-1a00-5bcc86e4d05c
md"""
## Q8
"""

# ╔═╡ 85f5171a-8a5e-11eb-3a2d-b5a338582dc8
q8b = [1 1]'

# ╔═╡ 8eef27de-8a5e-11eb-3bdc-150b09eb3971
q8a1 = [1 0]'

# ╔═╡ 951d2b2e-8a5e-11eb-15e4-9f286f08ed5a
q8a2 = [1 2]'

# ╔═╡ 9a8f05fa-8a5e-11eb-0712-6ddf30b0c038
q8p1 = p(q8a1, q8b)

# ╔═╡ c1aec63e-8a5e-11eb-3787-bbdecdf1fde1
q8p2 = p(q8a2, q8b)

# ╔═╡ bd57c16c-8a5e-11eb-0420-dfb1d3cd7f18
q8p1 + q8p2 # not equal to b

# ╔═╡ ea26e9a2-8a5e-11eb-107a-df619746d260
md"""
## Q9
"""

# ╔═╡ f90cf5e2-8a5e-11eb-39d2-2f93228b99e5
q9A = hcat(q8a1, q8a2)

# ╔═╡ 0663973a-8a5f-11eb-2300-c33617bb06b1
q9P = P(q9A)

# ╔═╡ 13868a14-8a5f-11eb-3db9-817cb263be30
md"""
## Q10
"""

# ╔═╡ 32184e7c-8a5f-11eb-127f-93f6e56481e8
q10a1 = [1 0]'

# ╔═╡ 37c8e3ea-8a5f-11eb-37f3-59e8ae99ff5f
q10a2 = [1 2]'

# ╔═╡ f3e17340-8a5e-11eb-08ed-f942fbaacd20
q10P1 = P(q10a1)

# ╔═╡ 5cb59162-8a5f-11eb-1ed6-373891fc90aa
q10P2 = P(q10a2)

# ╔═╡ 630d0932-8a5f-11eb-2803-dd7a28b649bf
q10P1P2 = q10P1 * q10P2

# ╔═╡ 76eacabe-8a5f-11eb-0994-2795f11b9068
q10P1P2^2 == q10P1P2

# ╔═╡ 6de08b90-8a5f-11eb-15bb-290d8c932fcf
q10P1P2* q10a1

# ╔═╡ f04303ec-8a5f-11eb-1693-b5e16cd3d677
md"""
## Q11
"""

# ╔═╡ 0bce6f98-8a60-11eb-00bb-7daf8b2c8315
q11A1 = [
	1 1;
	0 1;
	0 0
]

# ╔═╡ 72cc0156-8a60-11eb-1854-2546bb59784a
q11b1 = [2 3 4]'

# ╔═╡ 7949782c-8a60-11eb-06fb-2dac3de5fb58
q11A2 = [
	1 1;
	1 1;
	0 1
]

# ╔═╡ 861bbb46-8a60-11eb-3439-f5ea284a7058
q11b2 = [4 4 6]'

# ╔═╡ 9110e3e6-8a60-11eb-2793-5d883f950b1a
q11p1 = p(q11A1, q11b1)

# ╔═╡ a5dc48ee-8a60-11eb-1f99-f55ac517dd84
q11p2 = p(q11A2, q11b2)

# ╔═╡ bd12a134-8a60-11eb-15c2-bdc6898d88f0
q11e1 = e(q11p1, q11b1)

# ╔═╡ d05662b2-8a60-11eb-3547-2316a9d148d2
q11e1' * q11A1

# ╔═╡ fbe05000-8a60-11eb-0ed2-c34b788da62b
1+1

# ╔═╡ e49442d0-8a60-11eb-362f-5da5aa0f9a3e


# ╔═╡ Cell order:
# ╠═113d38f4-899a-11eb-27ef-bd783b2a1ba2
# ╠═dbeaa654-8a55-11eb-21f7-477091d9849b
# ╠═e71e2b92-8a56-11eb-3fe1-09ec3fe2771c
# ╠═f5f02be8-8a56-11eb-1228-35ab779e3be2
# ╠═0fb707e0-8a57-11eb-250a-959e36596fdb
# ╠═6a461b84-8a56-11eb-2a35-912dfb4ca40d
# ╠═7ec7c1fc-8a56-11eb-2866-c39711f9787e
# ╠═832c7ac6-8a56-11eb-245e-fdb85294c95b
# ╠═b6f7a2f4-8a56-11eb-0d45-29ecacfd5b93
# ╠═2296918c-8a57-11eb-25ca-7f1d9b3a6d6f
# ╠═93e57e9e-8a56-11eb-0dd8-cd087d26a743
# ╠═2bd4456e-8a57-11eb-0d27-87a4529d9b3c
# ╠═5f38e6bc-8a57-11eb-11ff-d101b514d5f5
# ╠═6e407d1e-8a57-11eb-3da2-0b5ca81c87a1
# ╠═754de3f8-8a57-11eb-0871-31b59419fe80
# ╟─9030867c-8a59-11eb-0dd4-233a23acba24
# ╟─8541822c-8a59-11eb-17eb-5936e9f7388d
# ╠═7d481752-8a57-11eb-33ed-dd651cc0f199
# ╠═8a2fa55c-8a57-11eb-2ac7-3d295325c893
# ╠═a3971ff4-8a57-11eb-2382-a1da4802eb2c
# ╠═17325ac8-8a58-11eb-3963-69200384e81f
# ╠═1ee9963a-8a58-11eb-338a-cf27c7fe7d74
# ╠═2fabe29a-8a58-11eb-0cb9-950bf424bda0
# ╠═4335141e-8a58-11eb-3f07-f74ccfd631d9
# ╠═9144618a-8a58-11eb-0664-bb1f527b3a93
# ╠═f5696b14-8a59-11eb-1185-57cecbd9e437
# ╠═33989bbc-8a5a-11eb-0d8e-39c297ea0627
# ╠═877bb944-8a5a-11eb-1cf8-414885322af4
# ╠═9281c784-8a5a-11eb-3cb7-179e43897555
# ╠═b13d1144-8a5a-11eb-113b-c15cd1ec4a65
# ╠═b97a79c6-8a5a-11eb-0009-a3b0a3dc6b45
# ╠═c3c2c33e-8a5a-11eb-2d7f-47784785ff89
# ╠═f69a69b0-8a5a-11eb-25d9-6374163cd317
# ╠═0d172a2a-8a5b-11eb-35bd-958247b6a94c
# ╠═18b99b60-8a5b-11eb-35c7-a30bf4517452
# ╠═2651f81c-8a5b-11eb-07b0-df1b53c0f16f
# ╠═5d4e180a-8a5b-11eb-3658-bb32f621f04b
# ╠═67ee890c-8a5b-11eb-2845-51dfba63507c
# ╠═827e5ad6-8a5b-11eb-13da-cdfddadb87da
# ╠═803ddbbe-8a5b-11eb-18da-b7da029a7fb0
# ╠═0eed69da-8a5c-11eb-3313-052ddd989a8f
# ╠═f3edfa3c-8a5b-11eb-23a1-0ff1b2a8c811
# ╠═073eb64e-8a5c-11eb-1609-4799913b89ff
# ╠═767e04de-8a5b-11eb-3645-2d2ce2984b8c
# ╠═44b4218a-8a5c-11eb-212f-9f128a8ec531
# ╠═5a9d18ce-8a5b-11eb-28c7-014c17246015
# ╠═6b5683b4-8a5c-11eb-1410-13b41f925244
# ╠═515a7384-8a5b-11eb-1ee9-8d10c0d2a66e
# ╠═8c2271a2-8a5c-11eb-3e47-53ce1261473e
# ╠═a1f04194-8a5c-11eb-33fc-a3ce4730ce7f
# ╠═b7d0ad6e-8a5c-11eb-0c90-0197079e8e9c
# ╠═c36aa0ee-8a5c-11eb-0806-d3f0366148eb
# ╠═ce18fc02-8a5c-11eb-004d-e11309767ac4
# ╠═dd96491e-8a5c-11eb-24a6-2f4c230c0311
# ╠═f45ac2aa-8a5d-11eb-2913-593cf0569cd3
# ╠═04281bc4-8a5e-11eb-2cd9-774a44c2a49f
# ╠═ed83a4b0-8a5d-11eb-00c1-d575a68861ba
# ╠═eac9248e-8a5d-11eb-35cb-8585b1265fb5
# ╠═4c449b28-8a5e-11eb-36f0-77625dc02ded
# ╠═7149ec0c-8a5e-11eb-1a00-5bcc86e4d05c
# ╠═85f5171a-8a5e-11eb-3a2d-b5a338582dc8
# ╠═8eef27de-8a5e-11eb-3bdc-150b09eb3971
# ╠═951d2b2e-8a5e-11eb-15e4-9f286f08ed5a
# ╠═9a8f05fa-8a5e-11eb-0712-6ddf30b0c038
# ╠═c1aec63e-8a5e-11eb-3787-bbdecdf1fde1
# ╠═bd57c16c-8a5e-11eb-0420-dfb1d3cd7f18
# ╠═ea26e9a2-8a5e-11eb-107a-df619746d260
# ╠═f90cf5e2-8a5e-11eb-39d2-2f93228b99e5
# ╠═0663973a-8a5f-11eb-2300-c33617bb06b1
# ╠═13868a14-8a5f-11eb-3db9-817cb263be30
# ╠═32184e7c-8a5f-11eb-127f-93f6e56481e8
# ╠═37c8e3ea-8a5f-11eb-37f3-59e8ae99ff5f
# ╠═f3e17340-8a5e-11eb-08ed-f942fbaacd20
# ╠═5cb59162-8a5f-11eb-1ed6-373891fc90aa
# ╠═630d0932-8a5f-11eb-2803-dd7a28b649bf
# ╠═76eacabe-8a5f-11eb-0994-2795f11b9068
# ╠═6de08b90-8a5f-11eb-15bb-290d8c932fcf
# ╠═f04303ec-8a5f-11eb-1693-b5e16cd3d677
# ╠═0bce6f98-8a60-11eb-00bb-7daf8b2c8315
# ╠═72cc0156-8a60-11eb-1854-2546bb59784a
# ╠═7949782c-8a60-11eb-06fb-2dac3de5fb58
# ╠═861bbb46-8a60-11eb-3439-f5ea284a7058
# ╠═9110e3e6-8a60-11eb-2793-5d883f950b1a
# ╠═a5dc48ee-8a60-11eb-1f99-f55ac517dd84
# ╠═bd12a134-8a60-11eb-15c2-bdc6898d88f0
# ╠═d05662b2-8a60-11eb-3547-2316a9d148d2
# ╠═fbe05000-8a60-11eb-0ed2-c34b788da62b
# ╠═e49442d0-8a60-11eb-362f-5da5aa0f9a3e
