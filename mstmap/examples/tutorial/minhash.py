import tmap as tm


def main():
    """ Main function """

    enc = tm.Minhash()

    mh_a = enc.from_binary_array(tm.VectorUchar([1, 1, 1, 1, 0, 1, 0, 1, 1, 0]))
    mh_b = enc.from_binary_array(tm.VectorUchar([1, 0, 1, 1, 0, 1, 1, 0, 1, 0]))
    mh_c = enc.from_binary_array(tm.VectorUchar([1, 0, 1, 1, 1, 1, 1, 0, 1, 0]))

    dist_a_b = enc.get_distance(mh_a, mh_b)
    dist_b_c = enc.get_distance(mh_b, mh_c)

    print(dist_a_b)
    print(dist_b_c)


if __name__ == "__main__":
    main()
