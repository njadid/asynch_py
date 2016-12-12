#ifndef ENDIAN_H
#define ENDIAN_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


typedef unsigned long long uint64_t;
typedef unsigned int uint32_t;
typedef unsigned short uint16_t;
typedef unsigned int uint8_t;


#if defined(BIG_ENDIAN) && !defined(LITTLE_ENDIAN)

#define htons(A) (A)
#define htonl(A) (A)
#define ntohs(A) (A)
#define ntohl(A) (A)

#define htobe16(x) (x)
#define htole16(x) ((((((uint16_t)(x)) >> 8))|((((uint16_t)(x)) << 8)))
#define be16toh(x) (x)
#define le16toh(x) ((((((uint16_t)(x)) >> 8))|((((uint16_t)(x)) << 8)))

#define htobe32(x) (x)
#define htole32(x) (((uint32_t)htole16(((uint16_t)(((uint32_t)(x)) >> 16)))) | (((uint32_t)htole16(((uint16_t)(x)))) << 16))
#define be32toh(x) (x)
#define le32toh(x) (((uint32_t)le16toh(((uint16_t)(((uint32_t)(x)) >> 16)))) | (((uint32_t)le16toh(((uint16_t)(x)))) << 16))

#define htobe64(x) (x)
#define htole64(x) (((uint64_t)htole32(((uint32_t)(((uint64_t)(x)) >> 32)))) | (((uint64_t)htole32(((uint32_t)(x)))) << 32))
#define be64toh(x) (x)
#define le64toh(x) (((uint64_t)le32toh(((uint32_t)(((uint64_t)(x)) >> 32)))) | (((uint64_t)le32toh(((uint32_t)(x)))) << 32))

#elif defined(LITTLE_ENDIAN) && !defined(BIG_ENDIAN)

#define htons(A) ((((uint16_t)(A) & 0xff00) >> 8) | \
(((uint16_t)(A) & 0x00ff) << 8))
#define htonl(A) ((((uint32_t)(A) & 0xff000000) >> 24) | \
(((uint32_t)(A) & 0x00ff0000) >> 8) | \
(((uint32_t)(A) & 0x0000ff00) << 8) | \
(((uint32_t)(A) & 0x000000ff) << 24))
#define ntohs htons
#define ntohl htonl

#define htobe16(x) htons(x)
#define htole16(x) (x)
#define be16toh(x) ntohs(x)
#define le16toh(x) (x)

#define htobe32(x) htonl(x)
#define htole32(x) (x)
#define be32toh(x) ntohl(x)
#define le32toh(x) (x)

#define htobe64(x) (((uint64_t)htonl(((uint32_t)(((uint64_t)(x)) >> 32)))) | (((uint64_t)htonl(((uint32_t)(x)))) << 32))
#define htole64(x) (x)
#define be64toh(x) (((uint64_t)ntohl(((uint32_t)(((uint64_t)(x)) >> 32)))) | (((uint64_t)ntohl(((uint32_t)(x)))) << 32))
#define le64toh(x) (x)

#else

#error "Either BIG_ENDIAN or LITTLE_ENDIAN must be #defined, but not both."

#endif

#endif //ENDIAN_H