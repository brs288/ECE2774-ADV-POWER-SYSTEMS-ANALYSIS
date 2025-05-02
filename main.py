import Validations

"""Verify test scenarios"""

# circ.add_solar("Solar_Scenario1", "bus7", 50, 0.95, "lagging")
Validations.Testing_Scenario1()

# circ.add_solar("Solar_Scenario2", "bus7", 20, 1, "unity") ... also sweep
Validations.Testing_Scenario2()

"""Verify edge and fail cases"""

# These can be set to True for testing purposes
fail_case_1 = False
fail_case_2 = False

circ1 = Validations.CreateSevenPowerBusSystem()
circ2 = Validations.CreateSevenPowerBusSystem()
circ3 = Validations.CreateSevenPowerBusSystem()

# Test 1: Solar generates no power
print("---EDGE CASE 1: 0 MW 0 MVAR---\n\n")
circ1.add_solar("Edge_Case1", "bus7", 0, 1, "unity")
Validations.NewtonRaphValidation(circ1)

# Test 2: Assign solar to slack
print("---EDGE CASE 2: 50 MW at bus1 ---\n\n")
circ2.add_solar("Edge_Case2", "bus1", 50, 1, "unity")
Validations.NewtonRaphValidation(circ2)

# Test 3: Add solar to bus that does not exist
print("---EDGE CASE 3: Add solar to non-existent bus---\n\n")
circ3.add_solar("Edge_Case3", "EDGE_CASE_3", 50, 1, "unity")

if fail_case_1:
    circ = Validations.CreateSevenPowerBusSystem()

    # Test 1: Make solar with power factor outside 0 < pf <= 1
    print("---FAIL CASE 1: 0 MW 50 MVAR---\n\n")
    circ.add_solar("Fail_Case1", "bus7", 50, -1, "lagging")
    Validations.NewtonRaphValidation(circ)

if fail_case_2:
    circ = Validations.CreateSevenPowerBusSystem()

    # Test 2: Faulty power factor mode
    print("---FAIL CASE 2: Faulty power factor mode---\n\n")
    circ.add_solar("Fail_Case2", "bus7", 50, 0.95, "FAIL CASE 2")
    Validations.NewtonRaphValidation(circ)